// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/sim/tracker_sim.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include "delphipp/reco/fwd_chambers.hpp"
#include "delphipp/reco/id.hpp"
#include "delphipp/reco/od.hpp"
#include "delphipp/reco/stw.hpp"
#include "delphipp/reco/tpc.hpp"
#include "delphipp/reco/vd.hpp"

namespace delphipp::sim {

namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;

float norm_phi(float phi) noexcept {
  while (phi < 0.0f) {
    phi += kTwoPi;
  }
  while (phi >= kTwoPi) {
    phi -= kTwoPi;
  }
  return phi;
}

float wrap_dphi(float a, float b) noexcept {
  float d = a - b;
  while (d > kPi) {
    d -= kTwoPi;
  }
  while (d < -kPi) {
    d += kTwoPi;
  }
  return d;
}

std::int32_t ascii_to_word(const char (&s)[5]) noexcept {
  return static_cast<std::int32_t>(
      (static_cast<std::uint32_t>(static_cast<unsigned char>(s[0])) << 24) |
      (static_cast<std::uint32_t>(static_cast<unsigned char>(s[1])) << 16) |
      (static_cast<std::uint32_t>(static_cast<unsigned char>(s[2])) << 8) |
      (static_cast<std::uint32_t>(static_cast<unsigned char>(s[3]))));
}

io::ParsedZebraBank make_bank(const char (&idh)[5],
                                std::int32_t orig_addr,
                                std::int32_t nl,
                                std::int32_t ns,
                                std::vector<std::int32_t> links,
                                std::vector<std::int32_t> words,
                                std::int32_t next_link = 0,
                                std::int32_t up_link = 0) {
  io::ParsedZebraBank b;
  b.idh = idh;
  b.idh_raw = ascii_to_word(idh);
  b.orig_addr = orig_addr;
  b.nl = nl;
  b.ns = ns;
  b.nd = static_cast<std::int32_t>(words.size());
  b.next_link = next_link;
  b.up_link = up_link;
  b.links = std::move(links);
  b.iwords = std::move(words);
  b.fwords.reserve(b.iwords.size());
  for (const std::int32_t w : b.iwords) {
    b.fwords.push_back(std::bit_cast<float>(static_cast<std::uint32_t>(w)));
  }
  return b;
}

std::vector<std::int32_t> make_blocklet_words(std::int32_t det_id,
                                              const std::vector<std::int32_t>& payload) {
  std::vector<std::int32_t> words;
  words.reserve(4u + payload.size());
  words.push_back(static_cast<std::int32_t>(4u + payload.size()));
  words.push_back(static_cast<std::int32_t>(static_cast<std::uint32_t>(det_id & 0xFF) << 8));
  words.push_back(0);
  words.push_back(0);
  words.insert(words.end(), payload.begin(), payload.end());
  return words;
}

std::uint32_t clamp_u32(long long val, long long lo, long long hi) noexcept {
  return static_cast<std::uint32_t>(std::clamp(val, lo, hi));
}

void encode_vd_hit(const SteppedHit& h, std::vector<std::int32_t>& out) {
  std::uint32_t layer_idx = 0;
  float best_dr = std::fabs(h.r_cm - reco::vd::kLayerRadiiCm[0]);
  for (std::uint32_t i = 1; i < 3u; ++i) {
    const float dr = std::fabs(h.r_cm - reco::vd::kLayerRadiiCm[i]);
    if (dr < best_dr) {
      best_dr = dr;
      layer_idx = i;
    }
  }
  const float r_layer = reco::vd::kLayerRadiiCm[layer_idx];
  const float sector_w = kTwoPi / 24.0f;
  const float phi_n = norm_phi(h.phi_rad);
  // Choose module so mod_phi <= phi_n (mod 2*pi) so strip * pitch >= 0
  const float shifted = norm_phi(phi_n - 0.5f * sector_w);
  const std::uint32_t mod_idx =
      static_cast<std::uint32_t>(std::floor(shifted / sector_w)) % 24u;
  const float mod_phi = (static_cast<float>(mod_idx) + 0.5f) * sector_w;
  const float dphi = norm_phi(phi_n - mod_phi);

  const reco::vd::VdConfig cfg{};
  const std::uint32_t strip = clamp_u32(
      std::llround(static_cast<double>(dphi * r_layer) / static_cast<double>(cfg.strip_pitch_cm)),
      0LL, 1023LL);
  const std::uint32_t adc = clamp_u32(
      std::llround(static_cast<double>(h.dedx_mev) * 500.0), 20LL, 255LL);

  const std::uint32_t w_rphi =
      (0u << 25) | (layer_idx << 23) | (mod_idx << 18) | (strip << 8) | adc;
  out.push_back(static_cast<std::int32_t>(w_rphi));

  // Double-sided z-strip readout only within non-negative representable local range
  if (h.z_cm >= 0.0f && h.z_cm <= 5.115f) {
    const std::uint32_t z_strip = clamp_u32(
        std::llround(static_cast<double>(h.z_cm) / static_cast<double>(cfg.strip_pitch_cm)),
        0LL, 1023LL);
    const std::uint32_t w_z =
        (1u << 25) | (layer_idx << 23) | (mod_idx << 18) | (z_strip << 8) | adc;
    out.push_back(static_cast<std::int32_t>(w_z));
  }
}

void encode_id_hit(const SteppedHit& h, std::vector<std::int32_t>& out) {
  const std::uint32_t wire_idx = clamp_u32(
      std::llround((static_cast<double>(h.r_cm) - 11.8) / (10.5 / 23.0)), 0LL, 23LL);
  const float r_wire = 11.8f + static_cast<float>(wire_idx) * (10.5f / 23.0f);
  const float sector_w = kTwoPi / static_cast<float>(reco::id::kNumSectors);
  const float phi_n = norm_phi(h.phi_rad);
  const std::uint32_t sec_idx =
      static_cast<std::uint32_t>(std::floor(phi_n / sector_w)) %
      static_cast<std::uint32_t>(reco::id::kNumSectors);
  const float sec_phi = (static_cast<float>(sec_idx) + 0.5f) * sector_w;
  const float dphi = wrap_dphi(phi_n, sec_phi);
  const std::uint32_t lr_bit = (dphi < 0.0f) ? 1u : 0u;
  const double drift_cm = static_cast<double>(std::fabs(dphi) * r_wire);
  const double drift_ns = drift_cm / 0.0040;
  const std::uint32_t tdc = clamp_u32(std::llround(drift_ns / 2.0), 1LL, 1023LL);

  const std::uint32_t w = (lr_bit << 20) | (sec_idx << 15) | (wire_idx << 10) | tdc;
  out.push_back(static_cast<std::int32_t>(w));
}

void encode_stw_hit(const SteppedHit& h, std::vector<std::int32_t>& out) {
  std::uint32_t layer_idx = 0;
  float best_dr = std::fabs(h.r_cm - reco::stw::kLayerRadiiCm[0]);
  for (std::uint32_t i = 1; i < static_cast<std::uint32_t>(reco::stw::kNumLayers); ++i) {
    const float dr = std::fabs(h.r_cm - reco::stw::kLayerRadiiCm[i]);
    if (dr < best_dr) {
      best_dr = dr;
      layer_idx = i;
    }
  }
  const float r_layer = reco::stw::kLayerRadiiCm[layer_idx];
  const float straw_w = kTwoPi / static_cast<float>(reco::stw::kNumStrawsPerLayer);
  const float phi_n = norm_phi(h.phi_rad);
  const float shifted = norm_phi(phi_n - 0.5f * straw_w);
  const std::uint32_t straw_idx =
      static_cast<std::uint32_t>(std::floor(shifted / straw_w)) %
      static_cast<std::uint32_t>(reco::stw::kNumStrawsPerLayer);
  const float phi_straw = (static_cast<float>(straw_idx) + 0.5f) * straw_w;
  const float dphi = norm_phi(phi_n - phi_straw);
  const double drift_ns = static_cast<double>(dphi * r_layer) / 0.005;
  const std::uint32_t tdc = clamp_u32(std::llround(drift_ns / 2.0), 1LL, 255LL);
  const std::uint32_t ph = clamp_u32(
      std::llround(static_cast<double>(h.dedx_mev) * 5000.0), 20LL, 255LL);

  const std::uint32_t w = (ph << 19) | (layer_idx << 16) | (straw_idx << 8) | tdc;
  out.push_back(static_cast<std::int32_t>(w));
}

void encode_tpc_hit(const SteppedHit& h, std::vector<std::int32_t>& out) {
  const std::uint32_t row_idx = clamp_u32(
      std::llround((static_cast<double>(h.r_cm) - 35.0) / 5.0), 0LL, 15LL);
  const float sector_w = kTwoPi / 6.0f;
  const std::uint32_t sec_idx =
      static_cast<std::uint32_t>(std::floor(norm_phi(h.phi_rad) / sector_w)) % 6u;
  const std::uint32_t endplate = (h.z_cm < 0.0f) ? 1u : 0u;

  const reco::tpc::TpcCalibration calib{};
  const double abs_z = std::min(static_cast<double>(std::fabs(h.z_cm)),
                                static_cast<double>(calib.endplate_z_cm));
  const double drift_dist = static_cast<double>(calib.endplate_z_cm) - abs_z;
  const double drift_us = drift_dist / static_cast<double>(calib.drift_velocity_cm_per_us);
  const double time_bin_d =
      static_cast<double>(calib.abst0_bins) + drift_us / static_cast<double>(calib.time_bin_us);
  const std::uint32_t time_bin = clamp_u32(std::llround(time_bin_d), 2LL, 511LL);
  const std::uint32_t charge = clamp_u32(
      std::llround(static_cast<double>(h.dedx_mev) * 2000.0), 20LL, 255LL);

  const std::uint32_t w =
      (endplate << 28) | (sec_idx << 24) | (row_idx << 20) | (time_bin << 8) | charge;
  out.push_back(static_cast<std::int32_t>(w));
}

void encode_od_hit(const SteppedHit& h, std::vector<std::int32_t>& out) {
  const std::uint32_t layer_idx = clamp_u32(
      std::llround((static_cast<double>(h.r_cm) - 197.0) / 2.25), 0LL, 4LL);
  const float plank_w = kTwoPi / static_cast<float>(reco::od::kNumPlanks);
  const std::uint32_t plank_idx =
      static_cast<std::uint32_t>(std::floor(norm_phi(h.phi_rad) / plank_w)) %
      static_cast<std::uint32_t>(reco::od::kNumPlanks);

  constexpr double kTMidNs = 64.0;
  const double dt_ns = static_cast<double>(h.z_cm) / 20.0;
  const std::uint32_t tdc_a = clamp_u32(std::llround((kTMidNs - dt_ns) / 0.5), 1LL, 255LL);
  const std::uint32_t tdc_c = clamp_u32(std::llround((kTMidNs + dt_ns) / 0.5), 1LL, 255LL);

  const std::uint32_t w =
      (plank_idx << 19) | (layer_idx << 16) | (tdc_c << 8) | tdc_a;
  out.push_back(static_cast<std::int32_t>(w));
}

void encode_fwd_hit(const SteppedHit& h, std::vector<std::int32_t>& out) {
  const double x = static_cast<double>(h.r_cm) * std::cos(static_cast<double>(h.phi_rad));
  const double y = static_cast<double>(h.r_cm) * std::sin(static_cast<double>(h.phi_rad));
  const std::uint32_t endcap_bit = (h.z_cm < 0.0f) ? 1u : 0u;
  constexpr double kPiD = std::numbers::pi_v<double>;
  constexpr std::array<double, 3> kAngles = {0.0, kPiD / 3.0, -kPiD / 3.0};

  for (std::uint32_t plane = 0; plane < 3u; ++plane) {
    const double u_cm = x * std::cos(kAngles[plane]) + y * std::sin(kAngles[plane]);
    const std::uint32_t u_raw = clamp_u32(std::llround(u_cm / 0.25 + 512.0), 1LL, 1023LL);
    const std::uint32_t w = (endcap_bit << 12) | (plane << 10) | u_raw;
    out.push_back(static_cast<std::int32_t>(w));
  }
}

}  // namespace

io::ParsedFzEvent digitize_trackers(const SimEvent& ev,
                                    const MagneticField& bfield,
                                    std::uint64_t seed) {
  std::vector<std::int32_t> vd_payload;
  std::vector<std::int32_t> id_payload;
  std::vector<std::int32_t> stw_payload;
  std::vector<std::int32_t> tpc_payload;
  std::vector<std::int32_t> od_payload;
  std::vector<std::int32_t> fca_payload;
  std::vector<std::int32_t> fcb_payload;

  for (std::size_t idx = 0; idx < ev.particles.size(); ++idx) {
    SimParticle p = ev.particles[idx];
    if (p.vertex_cm[0] == 0.0f && p.vertex_cm[1] == 0.0f && p.vertex_cm[2] == 0.0f) {
      p.vertex_cm = ev.primary_vertex_cm;
    }
    const std::uint64_t p_seed =
        seed ^ (static_cast<std::uint64_t>(idx + 1u) * 0x9E3779B97F4A7C15ULL);
    const auto hits = step_particle_in_bfield(p, bfield, p_seed);
    for (const auto& h : hits) {
      if (h.det_name == "VD") {
        encode_vd_hit(h, vd_payload);
      } else if (h.det_name == "ID") {
        encode_id_hit(h, id_payload);
      } else if (h.det_name == "STW") {
        encode_stw_hit(h, stw_payload);
      } else if (h.det_name == "TPC") {
        encode_tpc_hit(h, tpc_payload);
      } else if (h.det_name == "OD") {
        encode_od_hit(h, od_payload);
      } else if (h.det_name == "FCA") {
        encode_fwd_hit(h, fca_payload);
      } else if (h.det_name == "FCB") {
        encode_fwd_hit(h, fcb_payload);
      }
    }
  }

  io::ParsedFzEvent out;
  out.run_number = ev.run;
  out.event_number = ev.event;
  out.has_header_event_number = true;
  out.top_bank = "TOP ";

  constexpr std::int32_t kAddrRaw = 500;
  constexpr std::int32_t kAddrTop = 600;

  auto b_vd = make_bank("RAWD", 100, 0, 0, {},
                        make_blocklet_words(reco::vd::kDetectorId, vd_payload), 0, 110);
  auto p_vd = make_bank("RVD ", 110, 5, 5, {0, 0, 0, 0, 100}, {15, 0}, 0, kAddrRaw);

  auto b_id = make_bank("RAWD", 120, 0, 0, {},
                        make_blocklet_words(reco::id::kDetectorId, id_payload), 0, 130);
  auto p_id = make_bank("RID ", 130, 5, 5, {0, 0, 0, 0, 120}, {15, 0}, 0, kAddrRaw);

  auto b_tpc = make_bank("RAWD", 140, 0, 0, {},
                         make_blocklet_words(reco::tpc::kDetectorId, tpc_payload), 0, 150);
  auto p_tpc = make_bank("RTPC", 150, 5, 5, {0, 0, 0, 0, 140}, {15, 0}, 0, kAddrRaw);

  auto b_od = make_bank("RAWD", 160, 0, 0, {},
                        make_blocklet_words(reco::od::kDetectorId, od_payload), 0, 170);
  auto p_od = make_bank("ROD ", 170, 5, 5, {0, 0, 0, 0, 160}, {15, 0}, 0, kAddrRaw);

  auto b_fca = make_bank("RAWD", 180, 0, 0, {},
                         make_blocklet_words(reco::fwd::kDetectorIdFca, fca_payload), 0, 190);
  auto p_fca = make_bank("RFCA", 190, 5, 5, {0, 0, 0, 0, 180}, {15, 0}, 0, kAddrRaw);

  auto b_fcb = make_bank("RAWD", 200, 0, 0, {},
                         make_blocklet_words(reco::fwd::kDetectorIdFcb, fcb_payload), 0, 210);
  auto p_fcb = make_bank("RFCB", 210, 5, 5, {0, 0, 0, 0, 200}, {15, 0}, 0, kAddrRaw);

  auto b_stw = make_bank("STWR", 220, 0, 0, {},
                         make_blocklet_words(reco::stw::kDetectorId, stw_payload), 0, kAddrRaw);

  std::vector<std::int32_t> raw_links(22, 0);
  raw_links[6 - 1] = 110;   // RVD
  raw_links[7 - 1] = 130;   // RID
  raw_links[8 - 1] = 150;   // RTPC
  raw_links[10 - 1] = 170;  // ROD
  raw_links[15 - 1] = 190;  // RFCA
  raw_links[18 - 1] = 210;  // RFCB
  std::vector<std::int32_t> raw_words(15, 0);
  raw_words[0] = 15;
  auto raw_bank = make_bank("RAW ", kAddrRaw, 22, 22, std::move(raw_links),
                            std::move(raw_words), 0, kAddrTop);
  auto top_bank = make_bank("TOP ", kAddrTop, 7, 7, {kAddrRaw, 0, 0, 0, 0, 0, 0}, {1});

  out.banks = {
      std::move(b_vd),  std::move(p_vd),  std::move(b_id),  std::move(p_id),
      std::move(b_tpc), std::move(p_tpc), std::move(b_od),  std::move(p_od),
      std::move(b_fca), std::move(p_fca), std::move(b_fcb), std::move(p_fcb),
      std::move(b_stw), std::move(raw_bank), std::move(top_bank),
  };
  return out;
}

io::ParsedFzEvent build_sim_truth_record(const SimEvent& ev) {
  io::ParsedFzEvent out;
  out.run_number = ev.run;
  out.event_number = ev.event;
  out.has_header_event_number = true;
  out.top_bank = "SIM_";

  constexpr std::int32_t kAddrSim = 100;
  constexpr std::int32_t kAddrSp = 200;
  constexpr std::int32_t kAddrStBase = 300;

  const std::int32_t first_st = ev.particles.empty() ? 0 : kAddrStBase;
  std::vector<std::int32_t> sp_words(6, 0);
  sp_words[1] = std::bit_cast<std::int32_t>(ev.primary_vertex_cm[0]);
  sp_words[2] = std::bit_cast<std::int32_t>(ev.primary_vertex_cm[1]);
  sp_words[3] = std::bit_cast<std::int32_t>(ev.primary_vertex_cm[2]);

  auto sim_bank = make_bank("SIM_", kAddrSim, 2, 2, {kAddrSp, 0}, {ev.run, ev.event});
  auto sp_bank = make_bank("SP  ", kAddrSp, 2, 1, {first_st, 0}, std::move(sp_words), 0, kAddrSim);

  out.banks.push_back(std::move(sim_bank));
  out.banks.push_back(std::move(sp_bank));

  for (std::size_t i = 0; i < ev.particles.size(); ++i) {
    const auto& p = ev.particles[i];
    const std::int32_t addr = kAddrStBase + static_cast<std::int32_t>(i) * 10;
    const std::int32_t next_addr =
        (i + 1u < ev.particles.size()) ? (addr + 10) : 0;

    const float p_mag = std::sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
    const float e_val = (p.energy > 0.0f) ? p.energy : std::hypot(p_mag, p.mass);

    std::vector<std::int32_t> st_words(13, 0);
    st_words[0] = std::bit_cast<std::int32_t>(static_cast<float>(p.pdg));
    st_words[1] = std::bit_cast<std::int32_t>(p.px);
    st_words[2] = std::bit_cast<std::int32_t>(p.py);
    st_words[3] = std::bit_cast<std::int32_t>(p.pz);
    st_words[4] = std::bit_cast<std::int32_t>(e_val);
    st_words[5] = std::bit_cast<std::int32_t>(p_mag);
    st_words[6] = std::bit_cast<std::int32_t>(static_cast<float>(p.charge));
    st_words[7] = std::bit_cast<std::int32_t>(p.mass);

    out.banks.push_back(
        make_bank("ST  ", addr, 5, 1, {0, kAddrSp, 0, 0, 0},
                  std::move(st_words), next_addr, kAddrSp));
  }

  return out;
}

}  // namespace delphipp::sim
