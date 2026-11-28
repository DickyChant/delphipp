// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/sim/calo_pid_sim.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace delphipp::sim {

namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;
constexpr float kSpeedOfLightCmPerNs = 29.9792458f;

constexpr std::int32_t kDetBrich = 7;
constexpr std::int32_t kDetHpc = 9;
constexpr std::int32_t kDetTof = 11;
constexpr std::int32_t kDetHac = 13;
constexpr std::int32_t kDetMub = 14;
constexpr std::int32_t kDetEmf = 15;
constexpr std::int32_t kDetSat = 16;
constexpr std::int32_t kDetMus = 17;
constexpr std::int32_t kDetRichPartition = 17;
constexpr std::int32_t kDetLum = 21;
constexpr std::int32_t kDetHofSim = 22;
constexpr std::int32_t kDetMuf = 26;
constexpr std::int32_t kDetFrich = 27;
constexpr std::int32_t kDetStc = 32;

constexpr float kRichLiquidIndex = 1.273f;

std::int32_t ascii_to_word(std::string_view s) noexcept {
  if (s.size() < 4) {
    return 0;
  }
  const auto b0 = static_cast<std::uint32_t>(static_cast<unsigned char>(s[0]));
  const auto b1 = static_cast<std::uint32_t>(static_cast<unsigned char>(s[1]));
  const auto b2 = static_cast<std::uint32_t>(static_cast<unsigned char>(s[2]));
  const auto b3 = static_cast<std::uint32_t>(static_cast<unsigned char>(s[3]));
  return static_cast<std::int32_t>((b0 << 24) | (b1 << 16) | (b2 << 8) | b3);
}

float norm_phi(float phi) noexcept {
  while (phi < 0.0f) {
    phi += kTwoPi;
  }
  while (phi >= kTwoPi) {
    phi -= kTwoPi;
  }
  return phi;
}

float mass_from_pdg(std::int32_t pdg, bool is_charged) noexcept {
  const std::int32_t ap = (pdg < 0) ? -pdg : pdg;
  switch (ap) {
    case 11:   return 0.000511f;
    case 13:   return 0.105658f;
    case 22:   return 0.0f;
    case 211:  return 0.139570f;
    case 321:  return 0.493677f;
    case 2212: return 0.938272f;
    case 2112: return 0.939565f;
    case 130:
    case 310:  return 0.497611f;
    default:   return is_charged ? 0.139570f : 0.0f;
  }
}

struct ParticleKinematics {
  float px = 0.0f;
  float py = 0.0f;
  float pz = 0.0f;
  float pt = 0.0f;
  float pmag = 0.0f;
  float energy = 0.0f;
  float beta = 1.0f;
  float theta = 0.0f;
  float theta_beam = 0.0f;
  float phi = 0.0f;
  bool is_em = false;
  bool is_muon = false;
  bool is_hadron = false;
  bool is_charged = false;
};

bool compute_kinematics(const SimShowerParticle& p, ParticleKinematics& out) noexcept {
  if (!std::isfinite(p.px) || !std::isfinite(p.py) || !std::isfinite(p.pz) ||
      !std::isfinite(p.energy) || p.energy <= 0.0f) {
    return false;
  }
  const float pt = std::hypot(p.px, p.py);
  const float pmag = std::hypot(pt, p.pz);
  if (!std::isfinite(pmag) || pmag < 1e-6f) {
    return false;
  }

  const std::int32_t ap = (p.pdg < 0) ? -p.pdg : p.pdg;
  const bool is_em = (ap == 11 || ap == 22);
  const bool is_muon = (ap == 13);
  const bool is_charged = (p.charge != 0) || (ap == 11) || (ap == 13);
  const bool is_hadron = (!is_em && !is_muon);

  const float energy = (p.energy > pmag) ? p.energy : pmag;
  const float m = mass_from_pdg(p.pdg, is_charged);
  const float beta_mass = pmag / std::hypot(pmag, m);
  const float beta = (p.energy > pmag * 1.0001f)
                         ? std::clamp(pmag / energy, 0.0f, 1.0f)
                         : std::clamp(beta_mass, 0.0f, 1.0f);

  const float cos_th = std::clamp(p.pz / pmag, -1.0f, 1.0f);
  const float theta = std::acos(cos_th);
  const float theta_beam = std::min(theta, kPi - theta);
  const float phi = norm_phi(std::atan2(p.py, p.px));

  out = ParticleKinematics{
      .px = p.px,
      .py = p.py,
      .pz = p.pz,
      .pt = pt,
      .pmag = pmag,
      .energy = energy,
      .beta = beta,
      .theta = theta,
      .theta_beam = theta_beam,
      .phi = phi,
      .is_em = is_em,
      .is_muon = is_muon,
      .is_hadron = is_hadron,
      .is_charged = is_charged,
  };
  return true;
}

io::ParsedZebraBank make_zebra_bank(std::string_view idh,
                                    std::int32_t orig_addr,
                                    std::int32_t nl,
                                    std::int32_t ns,
                                    std::int32_t up_link,
                                    std::int32_t next_link,
                                    const std::vector<std::int32_t>& links,
                                    const std::vector<std::int32_t>& words) {
  io::ParsedZebraBank b;
  b.idh = std::string(idh);
  b.idh_raw = ascii_to_word(idh);
  b.nl = nl;
  b.ns = ns;
  b.nd = static_cast<std::int32_t>(words.size());
  b.orig_addr = orig_addr;
  b.next_link = next_link;
  b.up_link = up_link;
  b.links = links;
  b.iwords = words;
  b.fwords.reserve(words.size());
  for (const std::int32_t w : words) {
    b.fwords.push_back(std::bit_cast<float>(static_cast<std::uint32_t>(w)));
  }
  return b;
}

std::vector<std::int32_t> wrap_standard_blocklet(std::int32_t detector_id,
                                                 std::span<const std::int32_t> hits) {
  std::vector<std::int32_t> out;
  out.reserve(4 + hits.size());
  out.push_back(static_cast<std::int32_t>(4 + hits.size()));
  out.push_back(static_cast<std::int32_t>(static_cast<std::uint32_t>(detector_id) << 8));
  out.push_back(0);
  out.push_back(0);
  out.insert(out.end(), hits.begin(), hits.end());
  return out;
}

std::vector<std::int32_t> wrap_hac_blocklet(std::span<const std::int32_t> hits) {
  std::vector<std::int32_t> out;
  out.reserve(5 + hits.size());
  out.push_back(static_cast<std::int32_t>(5 + hits.size()));
  out.push_back(static_cast<std::int32_t>(static_cast<std::uint32_t>(kDetHac) << 8));
  out.push_back(0x5A5A);
  out.push_back(0);
  out.insert(out.end(), hits.begin(), hits.end());
  out.push_back(0);
  return out;
}

void digitize_hpc(const ParticleKinematics& k,
                  std::mt19937_64& rng,
                  std::vector<std::int32_t>& words) {
  if (!k.is_em || k.theta < 0.70f || k.theta > (kPi - 0.70f)) {
    return;
  }

  const std::int32_t phi_idx =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 24.0f))), 0, 23);
  const float pt_safe = std::max(k.pt, 1e-4f);
  const float z_impact = std::clamp(228.0f * (k.pz / pt_safe), -245.0f, 245.0f);

  std::int32_t best_ring = 2;
  float best_dist = 1e9f;
  for (std::int32_t r = 0; r < 6; ++r) {
    const float base_z = -210.0f + (static_cast<float>(r) + 0.5f) * 84.0f;
    const float mid_z = base_z + 11.0f;
    const float d = std::fabs(z_impact - mid_z);
    if (d < best_dist) {
      best_dist = d;
      best_ring = r;
    }
  }
  const float base_z = -210.0f + (static_cast<float>(best_ring) + 0.5f) * 84.0f;
  const std::int32_t time_bin = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::lround((z_impact - base_z) / 0.0432f)), 0, 511);
  const std::int32_t module_0 = best_ring * 24 + phi_idx;

  // Longitudinal Gamma shower profile over 9 sampling layers (~2 X0 each, R=208..260 cm)
  constexpr float kB = 0.5f;
  constexpr float kEcGev = 0.015f;
  const float t_max = std::max(1.5f, std::log(std::max(k.energy, 0.1f) / kEcGev) - 0.5f);
  const float a_shape = 1.0f + kB * t_max;

  std::array<float, 9> weights{};
  float wsum = 0.0f;
  for (std::size_t l = 0; l < 9; ++l) {
    const float t_mid = (static_cast<float>(l) + 0.5f) * 2.0f;
    const float w = std::pow(t_mid, a_shape - 1.0f) * std::exp(-kB * t_mid);
    weights[l] = w;
    wsum += w;
  }
  if (wsum <= 0.0f) {
    return;
  }

  std::normal_distribution<float> stoch(1.0f, 0.04f);
  for (std::size_t l = 0; l < 9; ++l) {
    const float fluct = std::clamp(stoch(rng), 0.7f, 1.3f);
    const float e_layer = k.energy * (weights[l] / wsum) * fluct;
    std::int32_t adc_rem = static_cast<std::int32_t>(std::lround(e_layer / 0.025f));
    while (adc_rem > 0) {
      const std::int32_t chunk = std::min<std::int32_t>(adc_rem, 255);
      const std::uint32_t w = (static_cast<std::uint32_t>(module_0) << 21) |
                              (static_cast<std::uint32_t>(l) << 17) |
                              (static_cast<std::uint32_t>(time_bin) << 8) |
                              static_cast<std::uint32_t>(chunk);
      words.push_back(static_cast<std::int32_t>(w));
      adc_rem -= chunk;
    }
  }
}

void digitize_emf(const ParticleKinematics& k,
                  std::mt19937_64& rng,
                  std::vector<std::int32_t>& words) {
  if (!k.is_em || k.theta_beam < 0.17f || k.theta_beam > 0.68f) {
    return;
  }

  const std::uint32_t endcap_bit = (k.pz >= 0.0f) ? 0u : 1u;
  const float abs_pz = std::max(std::fabs(k.pz), 1e-4f);
  const float x_cm = 284.0f * k.px / abs_pz;
  const float y_cm = 284.0f * k.py / abs_pz;
  const std::int32_t ix0 =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround(x_cm / 5.6f)), -30, 30);
  const std::int32_t iy0 =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround(y_cm / 5.6f)), -30, 30);

  std::array<std::array<float, 3>, 3> grid{};
  grid[1][1] = 0.76f * k.energy;
  grid[0][1] = 0.045f * k.energy;
  grid[2][1] = 0.045f * k.energy;
  grid[1][0] = 0.045f * k.energy;
  grid[1][2] = 0.045f * k.energy;
  grid[0][0] = 0.015f * k.energy;
  grid[0][2] = 0.015f * k.energy;
  grid[2][0] = 0.015f * k.energy;
  grid[2][2] = 0.015f * k.energy;

  constexpr float kMaxBlockEnergy = 30.0f;
  if (grid[1][1] > kMaxBlockEnergy) {
    const float excess = grid[1][1] - kMaxBlockEnergy;
    grid[1][1] = kMaxBlockEnergy;
    grid[0][1] += 0.25f * excess;
    grid[2][1] += 0.25f * excess;
    grid[1][0] += 0.25f * excess;
    grid[1][2] += 0.25f * excess;
  }

  std::normal_distribution<float> stoch(1.0f, 0.03f);
  for (std::int32_t dx = -1; dx <= 1; ++dx) {
    for (std::int32_t dy = -1; dy <= 1; ++dy) {
      const float e_blk =
          grid[static_cast<std::size_t>(dx + 1)][static_cast<std::size_t>(dy + 1)] *
          std::clamp(stoch(rng), 0.8f, 1.2f);
      const std::int32_t adc =
          std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround(e_blk * 134.0f)), 0, 4095);
      if (adc < 7) {
        continue;
      }
      const auto ux = static_cast<std::uint32_t>(ix0 + dx + 32);
      const auto uy = static_cast<std::uint32_t>(iy0 + dy + 32);
      const std::uint32_t w = (endcap_bit << 24) | (uy << 18) | (ux << 12) |
                              static_cast<std::uint32_t>(adc);
      words.push_back(static_cast<std::int32_t>(w));
    }
  }
}

void digitize_hac(const ParticleKinematics& k,
                  std::mt19937_64& rng,
                  std::vector<std::int32_t>& words) {
  if (k.is_em || k.theta < 0.25f || k.theta > (kPi - 0.25f)) {
    return;
  }

  const std::int32_t itheta = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor((k.theta - 0.25f) / (2.64f / 64.0f))), 0, 63);
  const std::int32_t iphi = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 96.0f))), 0, 95);

  std::normal_distribution<float> stoch(1.0f, 0.08f);
  const float e_tot = k.is_muon ? 1.60f : (k.energy * std::clamp(stoch(rng), 0.7f, 1.3f));
  const std::array<float, 4> layer_frac =
      k.is_muon ? std::array<float, 4>{0.25f, 0.25f, 0.25f, 0.25f}
                : std::array<float, 4>{0.28f, 0.36f, 0.23f, 0.13f};

  for (std::size_t lay = 0; lay < 4; ++lay) {
    const float e_lay = e_tot * layer_frac[lay];
    std::int32_t adc_rem = static_cast<std::int32_t>(std::lround(e_lay / 0.17f));
    if (adc_rem < 2) {
      continue;
    }
    while (adc_rem > 0) {
      const std::int32_t chunk = std::min<std::int32_t>(adc_rem, 255);
      const std::uint32_t w = (static_cast<std::uint32_t>(lay) << 24) |
                              (static_cast<std::uint32_t>(iphi) << 16) |
                              (static_cast<std::uint32_t>(itheta) << 8) |
                              static_cast<std::uint32_t>(chunk);
      words.push_back(static_cast<std::int32_t>(w));
      adc_rem -= chunk;
    }
  }
}

void digitize_stc(const ParticleKinematics& k, std::vector<std::int32_t>& words) {
  if (!k.is_em || k.theta_beam < 0.025f || k.theta_beam > 0.20f) {
    return;
  }
  const std::uint32_t arm_bit = (k.pz >= 0.0f) ? 0u : 1u;
  const float r_cm = 220.0f * std::tan(k.theta_beam);
  const std::int32_t ring_0 =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround((r_cm - 10.0f) / 3.5f)), 0, 9);
  const std::int32_t sector_0 = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 16.0f))), 0, 15);
  const std::int32_t adc =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround(k.energy / 0.05f)), 1, 4095);
  const std::uint32_t w = (arm_bit << 20) |
                          (static_cast<std::uint32_t>(sector_0) << 16) |
                          (static_cast<std::uint32_t>(ring_0) << 12) |
                          static_cast<std::uint32_t>(adc);
  words.push_back(static_cast<std::int32_t>(w));
}

void digitize_sat_lum(const ParticleKinematics& k,
                      std::vector<std::int32_t>& sat_words,
                      std::vector<std::int32_t>& lum_words) {
  if (!k.is_em || k.theta_beam < 0.004f || k.theta_beam > 0.18f) {
    return;
  }
  const std::uint32_t arm_bit = (k.pz >= 0.0f) ? 0u : 1u;
  const float r_cm = 232.0f * std::tan(k.theta_beam);
  const std::int32_t ring_0 =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround((r_cm - 13.0f) / 3.0f)), 0, 7);
  const std::int32_t sector_0 = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 24.0f))), 0, 23);
  const std::int32_t adc =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround(k.energy / 0.05f)), 1, 4095);
  const std::uint32_t w = (arm_bit << 20) |
                          (static_cast<std::uint32_t>(sector_0) << 15) |
                          (static_cast<std::uint32_t>(ring_0) << 12) |
                          static_cast<std::uint32_t>(adc);
  const auto word_i = static_cast<std::int32_t>(w);
  if (k.theta_beam >= 0.025f) {
    sat_words.push_back(word_i);
  }
  lum_words.push_back(word_i);
}

void digitize_tof(const ParticleKinematics& k, std::vector<std::int32_t>& words) {
  if (!k.is_charged || k.theta < 0.70f || k.theta > (kPi - 0.70f)) {
    return;
  }
  const float phi_deg = k.phi * (180.0f / kPi);
  std::int32_t s0 = static_cast<std::int32_t>(std::lround((82.5f - phi_deg) / 15.0f)) % 24;
  if (s0 < 0) {
    s0 += 24;
  }
  const std::int32_t sector = ((k.pz >= 0.0f) ? 1 : 25) + s0;
  const std::int32_t counter_num = (sector - 1) * 4 + 1;
  const std::int32_t idx = counter_num - 1;
  const std::int32_t slot = (idx < 16) ? 12 : (idx / 16);
  const std::int32_t chan = idx % 16;

  const float path_cm = 311.25f / std::max(std::sin(k.theta), 0.1f);
  const float time_ns = path_cm / (std::max(k.beta, 0.1f) * kSpeedOfLightCmPerNs);
  const std::int32_t tdc =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround(time_ns / 0.5f)), 1, 65534);

  const std::uint32_t w = (static_cast<std::uint32_t>(slot) << 27) |
                          (static_cast<std::uint32_t>(chan) << 17) |
                          static_cast<std::uint32_t>(tdc);
  words.push_back(static_cast<std::int32_t>(w));
}

void digitize_hof(const ParticleKinematics& k,
                  std::array<std::array<std::uint32_t, 4>, 2>& quad_masks) {
  const bool muon_hit = k.is_muon && (k.theta_beam >= 0.15f && k.theta_beam <= 0.85f);
  const bool punch_hit =
      k.is_charged && (k.energy > 20.0f) && (k.theta_beam >= 0.18f && k.theta_beam <= 0.65f);
  if (!muon_hit && !punch_hit) {
    return;
  }
  const std::size_t cap_idx = (k.pz >= 0.0f) ? 0u : 1u;
  const std::size_t quad_idx =
      static_cast<std::size_t>(std::clamp<std::int32_t>(
          static_cast<std::int32_t>(std::floor(k.phi / (0.5f * kPi))), 0, 3));
  const float r_hof = 450.0f * std::tan(k.theta_beam);
  const std::int32_t counter_0 =
      std::clamp<std::int32_t>(static_cast<std::int32_t>((r_hof - 40.0f) / 14.0f), 0, 27);
  quad_masks[cap_idx][quad_idx] |= (1u << static_cast<std::uint32_t>(counter_0));
}

void digitize_mub(const ParticleKinematics& k, std::vector<std::int32_t>& words) {
  if (!k.is_muon || k.theta < 0.65f || k.theta > (kPi - 0.65f)) {
    return;
  }
  const std::int32_t phi_idx = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 24.0f))), 0, 23);
  const std::int32_t sector_0 = ((k.pz >= 0.0f) ? 0 : 24) + phi_idx;
  constexpr std::uint32_t kPos0 = 1u;
  constexpr std::array<std::uint32_t, 4> kLayers = {0u, 1u, 3u, 4u};

  for (const std::uint32_t lay0 : kLayers) {
    const std::uint32_t ltd_raw = 16u + lay0 * 4u;
    const std::uint32_t w = (static_cast<std::uint32_t>(sector_0) << 22) |
                            (lay0 << 19) | (kPos0 << 16) | ltd_raw;
    words.push_back(static_cast<std::int32_t>(w));
  }
}

void digitize_muf(const ParticleKinematics& k, std::vector<std::int32_t>& words) {
  if (!k.is_muon || k.theta_beam < 0.15f || k.theta_beam > 0.88f) {
    return;
  }
  const std::int32_t quad_in_cap = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 8.0f))), 0, 7);
  const std::int32_t quad_0 = ((k.pz >= 0.0f) ? 0 : 8) + quad_in_cap;
  const float r_muf = 450.0f * std::tan(k.theta_beam);
  const std::int32_t chamber_0 =
      std::clamp<std::int32_t>(static_cast<std::int32_t>(std::floor(r_muf / 18.0f)), 0, 21);
  constexpr std::uint32_t kTaRaw = 40u;

  for (std::uint32_t lay0 = 0; lay0 < 2; ++lay0) {
    const std::uint32_t w = (lay0 << 19) |
                            (static_cast<std::uint32_t>(chamber_0) << 14) |
                            (static_cast<std::uint32_t>(quad_0) << 10) | kTaRaw;
    words.push_back(static_cast<std::int32_t>(w));
  }
}

void digitize_mus(const ParticleKinematics& k, std::vector<std::int32_t>& words) {
  if (!k.is_muon || k.theta_beam < 0.45f || k.theta_beam > 1.35f) {
    return;
  }
  const std::int32_t phi_oct = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 8.0f))), 0, 7);
  const std::int32_t det_0 = ((k.pz >= 0.0f) ? 0 : 8) + phi_oct;
  constexpr std::uint32_t kElem0 = 20u;
  constexpr std::uint32_t kAdc = 80u;

  // Emit strip+tube pairs in plane 1 (layers 0 & 1) and plane 2 (layers 3 & 2)
  constexpr std::array<std::uint32_t, 4> kLayers = {0u, 1u, 3u, 2u};
  for (const std::uint32_t lay0 : kLayers) {
    const std::uint32_t w = (static_cast<std::uint32_t>(det_0) << 16) |
                            (lay0 << 14) | (kElem0 << 8) | kAdc;
    words.push_back(static_cast<std::int32_t>(w));
  }
}

void digitize_rich(const ParticleKinematics& k,
                   std::mt19937_64& rng,
                   std::vector<std::int32_t>& brich_words,
                   std::vector<std::int32_t>& frich_words) {
  if (!k.is_charged || k.beta <= 1e-4f) {
    return;
  }
  const float cos_tc = 1.0f / (kRichLiquidIndex * k.beta);
  if (cos_tc >= 1.0f || cos_tc <= -1.0f) {
    return;  // Below Cherenkov threshold
  }
  const float theta_c = std::acos(std::clamp(cos_tc, -1.0f, 1.0f));
  const float sin_tc = std::sin(theta_c);
  const std::int32_t n_photons = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::lround(18.0f * sin_tc * sin_tc)), 2, 14);
  const std::int32_t sector_0 = std::clamp<std::int32_t>(
      static_cast<std::int32_t>(std::floor(k.phi / (kTwoPi / 24.0f))), 0, 23);

  std::normal_distribution<float> tc_res(0.0f, 0.004f);

  if (k.theta >= 0.68f && k.theta <= (kPi - 0.68f)) {
    const float pt_safe = std::max(k.pt, 1e-4f);
    const float z_cm = std::clamp(130.0f * (k.pz / pt_safe), -128.0f, 127.0f);
    const auto z_byte = static_cast<std::uint32_t>(
        std::clamp<std::int32_t>(static_cast<std::int32_t>(std::lround(z_cm + 128.0f)), 0, 255));
    for (std::int32_t i = 0; i < n_photons; ++i) {
      const float tc_meas = theta_c + tc_res(rng);
      const std::int32_t tdc = std::clamp<std::int32_t>(
          static_cast<std::int32_t>(std::lround((tc_meas - 0.02f) / 0.001f)), 1, 639);
      const std::uint32_t w = (z_byte << 16) |
                              (static_cast<std::uint32_t>(sector_0) << 10) |
                              static_cast<std::uint32_t>(tdc);
      brich_words.push_back(static_cast<std::int32_t>(w));
    }
  }

  if (k.theta_beam >= 0.15f && k.theta_beam <= 0.72f) {
    const std::uint32_t z_sign = (k.pz >= 0.0f) ? 0u : 1u;
    for (std::int32_t i = 0; i < n_photons; ++i) {
      const float tc_meas = theta_c + tc_res(rng);
      const std::int32_t tdc = std::clamp<std::int32_t>(
          static_cast<std::int32_t>(std::lround((tc_meas - 0.02f) / 0.001f)), 1, 639);
      const std::uint32_t w = (z_sign << 15) |
                              (static_cast<std::uint32_t>(sector_0) << 10) |
                              static_cast<std::uint32_t>(tdc);
      frich_words.push_back(static_cast<std::int32_t>(w));
    }
  }
}

}  // namespace

io::ParsedFzEvent digitize_calo_muon_rich(std::int32_t run,
                                          std::int32_t event,
                                          std::span<const SimShowerParticle> particles,
                                          std::uint64_t seed) {
  std::mt19937_64 rng(seed);

  std::vector<std::int32_t> hpc_hits;
  std::vector<std::int32_t> emf_hits;
  std::vector<std::int32_t> hac_hits;
  std::vector<std::int32_t> stc_hits;
  std::vector<std::int32_t> sat_hits;
  std::vector<std::int32_t> lum_hits;
  std::vector<std::int32_t> tof_hits;
  std::array<std::array<std::uint32_t, 4>, 2> hof_masks{};
  std::vector<std::int32_t> mub_hits;
  std::vector<std::int32_t> muf_hits;
  std::vector<std::int32_t> mus_hits;
  std::vector<std::int32_t> brich_hits;
  std::vector<std::int32_t> frich_hits;

  for (const auto& p : particles) {
    ParticleKinematics kin{};
    if (!compute_kinematics(p, kin)) {
      continue;
    }
    digitize_hpc(kin, rng, hpc_hits);
    digitize_emf(kin, rng, emf_hits);
    digitize_hac(kin, rng, hac_hits);
    digitize_stc(kin, stc_hits);
    digitize_sat_lum(kin, sat_hits, lum_hits);
    digitize_tof(kin, tof_hits);
    digitize_hof(kin, hof_masks);
    digitize_mub(kin, mub_hits);
    digitize_muf(kin, muf_hits);
    digitize_mus(kin, mus_hits);
    digitize_rich(kin, rng, brich_hits, frich_hits);
  }

  io::ParsedFzEvent ev;
  ev.run_number = run;
  ev.event_number = event;
  ev.has_header_event_number = true;
  ev.top_bank = "TOP ";

  std::int32_t next_addr = 100;
  auto alloc_addr = [&next_addr]() noexcept {
    const std::int32_t addr = next_addr;
    next_addr += 10;
    return addr;
  };

  const std::int32_t raw_addr = alloc_addr();
  const std::int32_t top_addr = alloc_addr();
  std::vector<std::int32_t> raw_links(22, 0);

  auto add_partition = [&](std::size_t raw_link_1based,
                           std::string_view pmb_name,
                           std::int32_t det_id,
                           const std::vector<std::int32_t>& hits,
                           bool is_hac) {
    if (hits.empty()) {
      return;
    }
    const std::int32_t rawd_addr = alloc_addr();
    const std::int32_t pmb_addr = alloc_addr();
    const auto blocklet = is_hac ? wrap_hac_blocklet(hits)
                                 : wrap_standard_blocklet(det_id, hits);
    ev.banks.push_back(
        make_zebra_bank("RAWD", rawd_addr, 0, 0, pmb_addr, 0, {}, blocklet));
    ev.banks.push_back(make_zebra_bank(
        pmb_name, pmb_addr, 5, 5, raw_addr, 0, {0, 0, 0, 0, rawd_addr}, {15, 0}));
    raw_links[raw_link_1based - 1] = pmb_addr;
  };

  add_partition(9, "RRIB", kDetRichPartition, brich_hits, false);
  add_partition(11, "RHPC", kDetHpc, hpc_hits, false);
  add_partition(12, "RTOF", kDetTof, tof_hits, false);
  add_partition(13, "RHAC", kDetHac, hac_hits, true);
  add_partition(14, "RMUB", kDetMub, mub_hits, false);
  add_partition(16, "RRIF", kDetRichPartition, frich_hits, false);
  add_partition(17, "RSAT", kDetSat, sat_hits, false);
  add_partition(19, "REMF", kDetEmf, emf_hits, false);
  add_partition(21, "RMUF", kDetMuf, muf_hits, false);
  add_partition(22, "RLUM", kDetLum, lum_hits, false);

  // HOF: one RAWD per hit endcap (module = kModuleSim = 22), chained via next_link
  std::vector<std::int32_t> hof_rawd_addrs;
  std::vector<std::vector<std::int32_t>> hof_rawd_words;
  for (std::size_t cap = 0; cap < 2; ++cap) {
    const auto& m = hof_masks[cap];
    if ((m[0] | m[1] | m[2] | m[3]) == 0u) {
      continue;
    }
    const std::int32_t numcap = static_cast<std::int32_t>(cap + 1);
    const std::uint32_t hdr2 =
        (static_cast<std::uint32_t>(numcap) << 16) |
        (static_cast<std::uint32_t>(kDetHofSim) << 8);
    hof_rawd_addrs.push_back(alloc_addr());
    hof_rawd_words.push_back({
        8,
        static_cast<std::int32_t>(hdr2),
        0,
        0,
        static_cast<std::int32_t>(m[0]),
        static_cast<std::int32_t>(m[1]),
        static_cast<std::int32_t>(m[2]),
        static_cast<std::int32_t>(m[3]),
    });
  }
  if (!hof_rawd_addrs.empty()) {
    const std::int32_t pmb_hof_addr = alloc_addr();
    for (std::size_t i = 0; i < hof_rawd_addrs.size(); ++i) {
      const std::int32_t next_hof =
          (i + 1 < hof_rawd_addrs.size()) ? hof_rawd_addrs[i + 1] : 0;
      ev.banks.push_back(make_zebra_bank(
          "RAWD", hof_rawd_addrs[i], 0, 0, pmb_hof_addr, next_hof, {}, hof_rawd_words[i]));
    }
    ev.banks.push_back(make_zebra_bank(
        "RHOF", pmb_hof_addr, 5, 5, raw_addr, 0, {0, 0, 0, 0, hof_rawd_addrs[0]}, {15, 0}));
    raw_links[20 - 1] = pmb_hof_addr;
  }

  // Legacy EDRA for MUS surround muon chambers (strip+tube pairs)
  if (!mus_hits.empty()) {
    ev.banks.push_back(make_zebra_bank(
        "EDRA", alloc_addr(), 0, 0, 0, 0, {}, wrap_standard_blocklet(kDetMus, mus_hits)));
  }

  // Legacy STRA for STIC forward luminosity calorimeter
  if (!stc_hits.empty()) {
    ev.banks.push_back(make_zebra_bank(
        "STRA", alloc_addr(), 0, 0, 0, 0, {}, wrap_standard_blocklet(kDetStc, stc_hits)));
  }

  ev.banks.push_back(make_zebra_bank(
      "RAW ", raw_addr, 22, 22, top_addr, 0, raw_links, {15, 0}));
  ev.banks.push_back(make_zebra_bank(
      "TOP ", top_addr, 7, 7, 0, 0, {raw_addr, 0, 0, 0, 0, 0, 0}, {1}));

  return ev;
}

}  // namespace delphipp::sim
