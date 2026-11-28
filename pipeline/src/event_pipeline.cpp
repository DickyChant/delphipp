// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/pipeline/event_pipeline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <unordered_set>
#include <vector>

#include "delphipp/helix.hpp"
#include "delphipp/io/raw_blocklet.hpp"
#include "delphipp/reco/ambiguity_solver.hpp"
#include "delphipp/reco/emf.hpp"
#include "delphipp/reco/fwd_chambers.hpp"
#include "delphipp/reco/had.hpp"
#include "delphipp/reco/hof.hpp"
#include "delphipp/reco/hpc.hpp"
#include "delphipp/reco/id.hpp"
#include "delphipp/reco/kalman_fitter.hpp"
#include "delphipp/reco/mub.hpp"
#include "delphipp/reco/muf.hpp"
#include "delphipp/reco/mus.hpp"
#include "delphipp/reco/od.hpp"
#include "delphipp/reco/rich.hpp"
#include "delphipp/reco/sat.hpp"
#include "delphipp/reco/stc.hpp"
#include "delphipp/reco/stw.hpp"
#include "delphipp/reco/tof.hpp"
#include "delphipp/reco/tpc.hpp"
#include "delphipp/reco/vd.hpp"
#include "delphipp/reco/vertex_fitter.hpp"
#include "delphipp/zebra_store.hpp"

namespace delphipp::pipeline {

namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;

float wrap_dphi(float a, float b) noexcept {
  float d = a - b;
  while (d > kPi) d -= kTwoPi;
  while (d < -kPi) d += kTwoPi;
  return d;
}

float norm_phi(float phi) noexcept {
  while (phi < 0.0f) phi += kTwoPi;
  while (phi >= kTwoPi) phi -= kTwoPi;
  return phi;
}

struct RawTrackHit {
  reco::tracking::CylinderMeasurement meas{};
  std::int32_t ter_id = -1;
  float dedx_sample = 0.0f;
};

struct RawCaloCluster {
  Vec3 xyz_cm{};
  float phi_rad = 0.0f;
  float theta_rad = 1.5707963f;
  float energy_gev = 0.0f;
  bool is_ecal = true;
  bool matched = false;
};

struct RawMuonHit {
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
};

struct RawRichMeasurement {
  std::int32_t detector_id = reco::rich::kDetectorIdBrich;
  float phi_rad = 0.0f;
  float theta_c_rad = 0.0f;
  std::int32_t num_photons = 0;
};

reco::muf::MufCalibration make_default_muf_calib() noexcept {
  reco::muf::MufCalibration cal{};
  cal.vanode_cm_per_ns = 25.0f;
  cal.wire_pitch_cm = 18.0f;
  for (std::size_t q = 0; q < reco::muf::kNumQuadrants; ++q) {
    const float angle = static_cast<float>(q % 4) * (0.5f * kPi);
    cal.qgeom[q][2] = 60.0f * std::cos(angle);
    cal.qgeom[q][3] = 60.0f * std::sin(angle);
    cal.svdrif[q] = {0.005f, 0.005f};
    cal.sx0par[q] = {0.5f, 0.5f};
  }
  return cal;
}

reco::mus::MusGeometry make_default_mus_geom() noexcept {
  reco::mus::MusGeometry geom{};
  for (std::size_t d = 0; d < reco::mus::kNumDetectors; ++d) {
    const float phi = (static_cast<float>(d % 8) + 0.5f) * (kTwoPi / 8.0f);
    const float z_off = (d < 8) ? 180.0f : -180.0f;
    for (std::size_t l = 0; l < reco::mus::kNumLayers; ++l) {
      geom.origin_ts1[d][l] = {-100.0f, -100.0f, 0.0f};
      geom.pitch_ts1[d][l] = {5.0f, 5.0f};
    }
    for (std::size_t p = 0; p < 2; ++p) {
      geom.local_to_delphi[d][p] = {
          470.0f * std::cos(phi), 470.0f * std::sin(phi), z_off,
          1.0f, 0.0f, 0.0f,
          0.0f, 1.0f, 0.0f,
          0.0f, 0.0f, 1.0f,
      };
    }
  }
  return geom;
}

PerigeeTrack track_from_p3_charge(float px,
                                  float py,
                                  float pz,
                                  std::int32_t charge,
                                  float bz,
                                  float eps = 0.0f,
                                  float z0 = 0.0f) {
  const float pt = std::hypot(px, py);
  const float p = std::hypot(pt, pz);
  const float theta =
      (p > 1e-6f) ? std::acos(std::clamp(pz / p, -1.0f, 1.0f)) : 1.5707963f;
  const float phi0 = std::atan2(py, px);
  const float signed_q = (charge == 0) ? 1.0f : static_cast<float>(charge);
  const float inv_r = (pt > 1e-5f) ? (signed_q * bz / pt) : 0.0f;
  return PerigeeTrack{
      .eps = eps,
      .z0 = z0,
      .theta = theta,
      .phi0 = phi0,
      .inv_r = inv_r,
  };
}

struct PartitionSpec {
  std::int32_t raw_link_idx;  // 1-based LQ(LRAW - il), per delsim36.car:23986
  const char* pmb_bank;
  const char* short_name;
  std::int32_t detector_id;
};

constexpr std::array<PartitionSpec, 18> kDelsimPartitions = {{
    {5,  "RCP ", "CP",   0},
    {6,  "RVD ", "VD",   3},
    {7,  "RID ", "ID",   4},
    {8,  "RTPC", "TPC",  5},
    {9,  "RRIB", "RIB",  17},
    {10, "ROD ", "OD",   8},
    {11, "RHPC", "HPC",  9},
    {12, "RTOF", "TOF",  11},
    {13, "RHAC", "HAC",  13},
    {14, "RMUB", "MUB",  14},
    {15, "RFCA", "FCA",  20},
    {16, "RRIF", "RIF",  17},
    {17, "RSAT", "SAT",  16},
    {18, "RFCB", "FCB",  25},
    {19, "REMF", "EMF",  15},
    {20, "RHOF", "HOF",  29},
    {21, "RMUF", "MUF",  26},
    {22, "RLUM", "LUM",  21},
}};

constexpr std::array<const char*, 18> kLegacyRawBankIds = {
    "EGRA", "EQRA", "EDRA", "EMRA", "EFRA", "RAW ",
    "VDRA", "IDRA", "TPRA", "ODRA", "HPRA", "HARA",
    "FCRA", "RIRA", "SARA", "STRA", "TORA", "MURA",
};

}  // namespace

std::string detector_name(std::int32_t detector_id) {
  switch (detector_id) {
    case 3:  return "VD";
    case 4:  return "ID";
    case 5:  return "TPC";
    case 8:  return "OD";
    case 9:  return "HPC";
    case 11: return "TOF";
    case 13: return "HAC";
    case 14: return "MUB";
    case 15: return "EMF";
    case 16: return "SAT";
    case 17: return "RICH";
    case 20: return "FCA";
    case 21: return "STC";
    case 25: return "FCB";
    case 26: return "MUF";
    case 29: return "HOF";
    case 30: return "MUS";
    case 31: return "STW";
    default: return "DET" + std::to_string(detector_id);
  }
}

EventReconstructionPipeline::EventReconstructionPipeline(
    const MagneticField& bfield, float abst0_val) noexcept
    : bfield_(bfield), abst0_val_(abst0_val) {}

EventRecoOutput EventReconstructionPipeline::reconstruct_event(
    const io::ParsedFzEvent& ev,
    std::int32_t fallback_run,
    std::int32_t fallback_evt) const {
  EventRecoOutput out;
  out.run = (ev.run_number != 0) ? ev.run_number : fallback_run;
  out.event = ev.has_header_event_number
                  ? ev.event_number
                  : ((ev.event_number > 0) ? ev.event_number : fallback_evt);
  out.top_bank = ev.top_bank;
  out.zebra_banks = ev.banks.size();

  const auto b0 = bfield_.evaluate({0.0f, 0.0f, 0.0f});

  std::size_t zstore_words = 1024;
  for (const auto& b : ev.banks) {
    zstore_words += static_cast<std::size_t>(
        std::max<std::int32_t>(b.nl, 0) + std::max<std::int32_t>(b.nd, 0) + 32);
  }
  zstore_words = std::max<std::size_t>(zstore_words, 65536);
  OwnedZebraStore zstore(zstore_words);
  auto zview = zstore.view();
  const auto bank_heads = ev.populate_zebra_store(zview);

  std::vector<dst::NeutralInput> neutrals;
  std::vector<RawTrackHit> raw_track_hits;
  std::vector<RawCaloCluster> raw_calo_clusters;
  std::vector<RawMuonHit> raw_muon_hits;
  std::vector<RawRichMeasurement> raw_rich_hits;
  bool is_bhabha = false;
  std::unordered_set<std::int32_t> processed_rawd_links;

  reco::hof::HofCalibration hof_calib{};
  reco::hof::HofStats hof_stats{};
  hof_stats.begin_event();
  const std::int32_t eqra_head = bank_heads.contains("EQRA") ? bank_heads.at("EQRA") : 0;

  reco::tpc::TpcCalibration tpc_calib{};
  if (abst0_val_ > 0.0f) {
    tpc_calib.abst0_bins = abst0_val_;
  }
  const auto muf_calib = make_default_muf_calib();
  const auto mus_geom = make_default_mus_geom();

  auto push_calo_xyz = [&](const Vec3& xyz, float energy_gev, bool is_ecal) {
    if (!std::isfinite(energy_gev) || energy_gev <= 0.0f) {
      return;
    }
    const float r = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
    if (r <= 1e-3f) {
      return;
    }
    raw_calo_clusters.push_back(RawCaloCluster{
        .xyz_cm = xyz,
        .phi_rad = norm_phi(std::atan2(xyz[1], xyz[0])),
        .theta_rad = std::acos(std::clamp(xyz[2] / r, -1.0f, 1.0f)),
        .energy_gev = energy_gev,
        .is_ecal = is_ecal,
        .matched = false,
    });
  };

  auto decode_partition_blocklet = [&](const std::string& pname,
                                       std::span<const std::int32_t> data) {
    if (pname == "VD") {
      const auto strips = reco::vd::unpack_raw_blocklet(data);
      const auto clusters = reco::vd::find_clusters(strips);
      const auto ter_ids = reco::vd::write_ter_banks(out.tanagra, clusters);
      for (std::size_t k = 0; k < clusters.size(); ++k) {
        const auto& c = clusters[k];
        const float r_cm = (c.layer >= 1 && c.layer <= 3)
                               ? reco::vd::kLayerRadiiCm[static_cast<std::size_t>(c.layer - 1)]
                               : reco::vd::kLayerRadiiCm[0];
        const float mod_phi = (static_cast<float>(c.module - 1) + 0.5f) * (kTwoPi / 24.0f);
        const float phi_rad = c.is_z_side ? mod_phi : norm_phi(mod_phi + c.local_coord_cm / r_cm);
        raw_track_hits.push_back(RawTrackHit{
            .meas = {
                .r_cm = r_cm,
                .phi_rad = phi_rad,
                .z_cm = c.is_z_side ? c.local_coord_cm : 0.0f,
                .sigma_rphi_cm = 0.002f,
                .sigma_z_cm = 0.005f,
                .has_z = c.is_z_side,
                .material_x0 = 0.004f,
            },
            .ter_id = static_cast<std::int32_t>(ter_ids[k]),
        });
      }
    } else if (pname == "ID") {
      const auto hits = reco::id::unpack_raw_blocklet(data);
      const auto segs = reco::id::reconstruct_segments(hits);
      const auto ter_ids = reco::id::write_ter_banks(out.tanagra, segs);
      for (std::size_t k = 0; k < segs.size(); ++k) {
        const auto& s = segs[k];
        raw_track_hits.push_back(RawTrackHit{
            .meas = {
                .r_cm = s.r_ref_cm,
                .phi_rad = norm_phi(s.phi_rad),
                .z_cm = 0.0f,
                .sigma_rphi_cm = 0.015f,
                .sigma_z_cm = 2.0f,
                .has_z = false,
                .material_x0 = 0.01f,
            },
            .ter_id = static_cast<std::int32_t>(ter_ids[k]),
        });
      }
    } else if (pname == "TPC") {
      const auto pad_hits = reco::tpc::unpack_raw_blocklet(data, tpc_calib);
      const auto pts = reco::tpc::reconstruct_space_points(pad_hits, tpc_calib);
      const auto ter_ids = reco::tpc::write_ter_banks(out.tanagra, pts);
      for (std::size_t k = 0; k < pts.size(); ++k) {
        const auto& p = pts[k];
        if (p.r_cm <= 1.0f) continue;
        raw_track_hits.push_back(RawTrackHit{
            .meas = {
                .r_cm = p.r_cm,
                .phi_rad = norm_phi(p.r_phi_cm / p.r_cm),
                .z_cm = p.z_cm,
                .sigma_rphi_cm = 0.025f,
                .sigma_z_cm = 0.09f,
                .has_z = true,
                .material_x0 = 0.003f,
            },
            .ter_id = static_cast<std::int32_t>(ter_ids[k]),
            .dedx_sample = p.charge * 0.025f,
        });
      }
    } else if (pname == "OD") {
      const auto tubes = reco::od::unpack_raw_blocklet(data);
      const auto elems = reco::od::reconstruct_planks(tubes);
      const auto ter_ids = reco::od::write_ter_banks(out.tanagra, elems);
      for (std::size_t k = 0; k < elems.size(); ++k) {
        const auto& e = elems[k];
        if (e.r_cm <= 1.0f) continue;
        raw_track_hits.push_back(RawTrackHit{
            .meas = {
                .r_cm = e.r_cm,
                .phi_rad = norm_phi(e.r_phi_cm / e.r_cm),
                .z_cm = e.z_cm,
                .sigma_rphi_cm = 0.015f,
                .sigma_z_cm = 3.5f,
                .has_z = true,
                .material_x0 = 0.015f,
            },
            .ter_id = static_cast<std::int32_t>(ter_ids[k]),
        });
      }
    } else if (pname == "FCA" || pname == "FCB") {
      const std::int32_t fdet = (pname == "FCB") ? reco::fwd::kDetectorIdFcb
                                               : reco::fwd::kDetectorIdFca;
      const auto wires = reco::fwd::unpack_raw_blocklet(data, fdet);
      for (std::int32_t endcap = 1; endcap <= 2; ++endcap) {
        std::vector<reco::fwd::ForwardWireHit> sub;
        for (const auto& w : wires) {
          if (w.endcap == endcap) sub.push_back(w);
        }
        if (const auto seg = reco::fwd::fit_segment(sub); seg.has_value()) {
          const auto ter_ids = reco::fwd::write_ter_banks(out.tanagra, std::span{&(*seg), 1}, fdet);
          const float r_xy = std::hypot(seg->x_cm, seg->y_cm);
          if (r_xy > 10.0f && r_xy < 260.0f) {
            raw_track_hits.push_back(RawTrackHit{
                .meas = {
                    .r_cm = r_xy,
                    .phi_rad = norm_phi(std::atan2(seg->y_cm, seg->x_cm)),
                    .z_cm = seg->z_cm,
                    .sigma_rphi_cm = 0.05f,
                    .sigma_z_cm = 0.5f,
                    .has_z = true,
                    .material_x0 = 0.02f,
                },
                .ter_id = static_cast<std::int32_t>(ter_ids.front()),
            });
          }
        }
      }
    } else if (pname == "HPC") {
      const auto hits = reco::hpc::unpack_raw_blocklet(data);
      const auto showers = reco::hpc::reconstruct_showers(hits);
      reco::hpc::write_ter_banks(out.tanagra, showers);
      for (const auto& s : showers) {
        push_calo_xyz({s.r_cm * std::cos(s.phi_rad), s.r_cm * std::sin(s.phi_rad), s.z_cm},
                      s.energy_gev, true);
      }
    } else if (pname == "EMF") {
      const auto blocks = reco::emf::unpack_raw_blocklet(data);
      const auto clusters = reco::emf::find_clusters(blocks);
      reco::emf::write_ter_banks(out.tanagra, clusters);
      for (const auto& c : clusters) {
        push_calo_xyz({c.x_cm, c.y_cm, c.z_cm}, c.energy_gev, true);
      }
    } else if (pname == "HAC") {
      const auto towers = reco::had::unpack_raw_blocklet(data);
      const auto clusters = reco::had::find_clusters(towers, 0.25f);
      reco::had::write_ter_banks(out.tanagra, clusters);
      for (const auto& cl : clusters) {
        const float st = std::sin(cl.theta_rad);
        push_calo_xyz({cl.r_cm * st * std::cos(cl.phi_rad),
                       cl.r_cm * st * std::sin(cl.phi_rad),
                       cl.r_cm * std::cos(cl.theta_rad)},
                      cl.energy_gev, false);
      }
    } else if (pname == "STC") {
      const auto towers = reco::stc::unpack_raw_blocklet(data);
      const auto clusters = reco::stc::find_clusters(towers);
      reco::stc::write_ter_banks(out.tanagra, clusters);
      is_bhabha = is_bhabha || reco::stc::is_bhabha_candidate(clusters);
      for (const auto& c : clusters) {
        push_calo_xyz({c.r_cm * std::cos(c.phi_rad), c.r_cm * std::sin(c.phi_rad), c.z_cm},
                      c.energy_gev, true);
      }
    } else if (pname == "SAT" || pname == "LUM") {
      const std::int32_t sdet = (pname == "LUM") ? reco::sat::kDetectorIdVsat
                                               : reco::sat::kDetectorIdSat;
      const auto hits = reco::sat::unpack_raw_blocklet(data);
      const auto clusters = reco::sat::reconstruct_clusters(hits);
      reco::sat::write_ter_banks(out.tanagra, clusters, sdet);
      for (const auto& c : clusters) {
        push_calo_xyz({c.r_cm * std::cos(c.phi_rad), c.r_cm * std::sin(c.phi_rad), c.z_cm},
                      c.energy_gev, true);
      }
    } else if (pname == "STW") {
      const auto hits = reco::stw::unpack_raw_blocklet(data);
      const auto elems = reco::stw::reconstruct_elements(hits);
      const auto ter_ids = reco::stw::write_ter_banks(out.tanagra, elems);
      for (std::size_t k = 0; k < elems.size(); ++k) {
        const auto& e = elems[k];
        raw_track_hits.push_back(RawTrackHit{
            .meas = {
                .r_cm = e.r_ref_cm,
                .phi_rad = norm_phi(e.phi_rad),
                .z_cm = 0.0f,
                .sigma_rphi_cm = 0.015f,
                .sigma_z_cm = 2.0f,
                .has_z = false,
                .material_x0 = 0.008f,
            },
            .ter_id = static_cast<std::int32_t>(ter_ids[k]),
        });
      }
    } else if (pname == "TOF") {
      const auto tof_hits = reco::tof::unpack_raw_blocklet(data);
      if (!tof_hits.empty()) {
        reco::tof::create_td_banks(out.tanagra, tof_hits);
      }
    } else if (pname == "MUB") {
      const auto hits = reco::mub::unpack_raw_blocklet(data);
      reco::mub::write_ter_banks(out.tanagra, hits);
      for (const auto& h : hits) {
        raw_muon_hits.push_back(RawMuonHit{.phi_rad = norm_phi(h.phi_rad), .z_cm = h.z_cm});
      }
    } else if (pname == "MUF") {
      const auto trips = reco::muf::unpack_raw_blocklet(data);
      const auto space_hits = reco::muf::reconstruct_hits(trips, muf_calib);
      reco::muf::write_ter_banks(out.tanagra, space_hits);
      for (const auto& h : space_hits) {
        const float z_cm = (h.quadrant <= 8) ? 450.0f : -450.0f;
        const float phi_rad = norm_phi(std::atan2(h.drift_dist_cm, h.x_wire_cm));
        raw_muon_hits.push_back(RawMuonHit{.phi_rad = phi_rad, .z_cm = z_cm});
      }
    } else if (pname == "MUS") {
      const auto raw_hits = reco::mus::unpack_raw_blocklet(data);
      const auto clusters = reco::mus::find_clusters(raw_hits);
      const auto sp = reco::mus::combine_clusters(clusters, mus_geom);
      const auto tes = reco::mus::build_te(sp, mus_geom);
      reco::mus::write_ter_banks(out.tanagra, tes);
      for (const auto& b : tes) {
        if (b.r > 1.0f) {
          raw_muon_hits.push_back(RawMuonHit{.phi_rad = norm_phi(b.r_phi / b.r), .z_cm = b.z});
        }
      }
    } else if (pname == "RIB" || pname == "RIF") {
      const std::int32_t rdet = (pname == "RIF") ? reco::rich::kDetectorIdFrich
                                               : reco::rich::kDetectorIdBrich;
      const auto photons = reco::rich::unpack_raw_blocklet(data, rdet);
      for (const auto& ph : photons) {
        raw_rich_hits.push_back(RawRichMeasurement{
            .detector_id = rdet,
            .phi_rad = norm_phi(ph.phi_rad),
            .theta_c_rad = ph.theta_c_rad,
            .num_photons = 1,
        });
      }
    }
  };

  auto process_rawd_chain = [&](const char* part_name,
                                std::int32_t default_det_id,
                                std::int32_t first_rawd) {
    if (first_rawd <= 0 || processed_rawd_links.contains(first_rawd)) {
      return;
    }
    const auto blocklets = io::read_raw_blocklet_chain(zview, first_rawd);
    for (const auto& bl : blocklets) {
      processed_rawd_links.insert(bl.bank_link);
    }
    if (blocklets.empty()) {
      return;
    }

    out.raw_blocklets += blocklets.size();
    out.detector_blocklet_counts[part_name] += blocklets.size();

    const std::string pname(part_name);
    if (pname == "HOF") {
      const auto hof_res = reco::hof::unpack_td(
          zview, first_rawd, eqra_head, hof_calib, hof_stats);
      out.hof_hits += hof_res.counter_hits;
      if (!hof_res.hits.empty()) {
        std::vector<float> flat;
        flat.reserve(hof_res.hits.size() * 8);
        for (const auto& h : hof_res.hits) {
          flat.insert(flat.end(), h.words.begin(), h.words.end());
          raw_muon_hits.push_back(RawMuonHit{
              .phi_rad = norm_phi(h.words[2] * 0.25f),
              .z_cm = (h.words[1] > 0.0f) ? 350.0f : -350.0f,
          });
        }
        out.tanagra.add_tdr(TdrBank{
            .detector_id = reco::hof::kDetectorId,
            .submod_id = 1,
            .descriptor = 290,
            .num_points = static_cast<std::int32_t>(hof_res.hits.size()),
            .words_per_point = 8,
            .data = std::move(flat),
        });
      }
    }

    for (const auto& bl : blocklets) {
      const std::int32_t det_id =
          (default_det_id > 0) ? default_det_id : bl.detector_id;
      const auto tdr_id = out.tanagra.add_tdr(TdrBank{
          .detector_id = det_id,
          .submod_id = bl.subdet_id,
          .descriptor = bl.flags,
          .num_points = 1,
          .words_per_point = static_cast<std::int32_t>(
              std::min<std::size_t>(bl.data.size(), 4)),
          .data = {
              static_cast<float>(bl.num_words),
              static_cast<float>(det_id),
              static_cast<float>(bl.subdet_id),
              abst0_val_,
          },
      });
      out.tanagra.add_ter(TerBank{
          .detector_id = det_id,
          .submod_id = bl.subdet_id,
          .descriptor = bl.flags,
          .iwords = {bl.num_words, det_id, bl.subdet_id},
          .fwords = {abst0_val_},
          .td_refs = {tdr_id},
      });
      decode_partition_blocklet(pname, bl.data);
    }
  };

  // 1. Traverse DELSIM TOP -> RAW -> 18 Partition Master Banks (LQ(LRAW-il)) -> RAWD (LQ(LPMB-5))
  std::int32_t lraw_root = 0;
  if (auto it = bank_heads.find("TOP "); it != bank_heads.end()) {
    const std::int32_t ltop = it->second;
    if (zview.bank_nl(ltop) >= 1) {
      lraw_root = zview.bank_link(ltop, 1);
    }
  }
  if (lraw_root <= 0) {
    if (auto it = bank_heads.find("RAW "); it != bank_heads.end()) {
      lraw_root = it->second;
    }
  }

  for (const auto& part : kDelsimPartitions) {
    std::int32_t lpmb = 0;
    if (lraw_root > 0 && zview.bank_nl(lraw_root) >= part.raw_link_idx) {
      lpmb = zview.bank_link(lraw_root, part.raw_link_idx);
    } else if (lraw_root <= 0) {
      if (auto it = bank_heads.find(part.pmb_bank); it != bank_heads.end()) {
        lpmb = it->second;
      }
    }
    while (lpmb > 0) {
      if (zview.bank_nl(lpmb) >= 5) {
        process_rawd_chain(
            part.short_name, part.detector_id, zview.bank_link(lpmb, 5));
      }
      lpmb = zview.bank_next(lpmb);
    }
  }

  // 2. Fallback: unpack HOF EGRA & legacy raw blocklet banks (flat *RA banks or unlinked RAWD)
  if (auto it = bank_heads.find("EGRA"); it != bank_heads.end()) {
    process_rawd_chain("HOF", reco::hof::kDetectorId, it->second);
  }
  auto part_from_legacy = [](std::string_view raw_id, std::int32_t det_id) -> std::string {
    if (raw_id == "EDRA") return "MUS";
    if (raw_id == "VDRA") return "VD";
    if (raw_id == "IDRA") return "ID";
    if (raw_id == "TPRA") return "TPC";
    if (raw_id == "ODRA") return "OD";
    if (raw_id == "HPRA") return "HPC";
    if (raw_id == "HARA") return "HAC";
    if (raw_id == "EMRA") return "EMF";
    if (raw_id == "FCRA") return "FCA";
    if (raw_id == "RIRA") return "RIB";
    if (raw_id == "SARA") return "SAT";
    if (raw_id == "STRA") return "STC";
    if (raw_id == "STWR") return "STW";
    if (raw_id == "TORA") return "TOF";
    if (raw_id == "MURA") return "MUB";
    return detector_name(det_id);
  };
  for (const char* raw_id : kLegacyRawBankIds) {
    if (std::string(raw_id) == "RAW ") {
      continue;
    }
    if (auto it = bank_heads.find(raw_id); it != bank_heads.end()) {
      const auto blocklets = io::read_raw_blocklet_chain(zview, it->second);
      for (const auto& bl : blocklets) {
        if (processed_rawd_links.contains(bl.bank_link)) {
          continue;
        }
        processed_rawd_links.insert(bl.bank_link);
        std::string dname = part_from_legacy(raw_id, bl.detector_id);
        if (std::string(raw_id) == "EGRA" &&
            bl.detector_id == reco::hof::kModuleSim) {
          dname = "HOF";
        }
        ++out.raw_blocklets;
        ++out.detector_blocklet_counts[dname];
        const auto tdr_id = out.tanagra.add_tdr(TdrBank{
            .detector_id = bl.detector_id,
            .submod_id = bl.subdet_id,
            .descriptor = bl.flags,
            .num_points = 1,
            .words_per_point = static_cast<std::int32_t>(
                std::min<std::size_t>(bl.data.size(), 4)),
            .data = {
                static_cast<float>(bl.num_words),
                static_cast<float>(bl.detector_id),
                static_cast<float>(bl.subdet_id),
                abst0_val_,
            },
        });
        out.tanagra.add_ter(TerBank{
            .detector_id = bl.detector_id,
            .submod_id = bl.subdet_id,
            .descriptor = bl.flags,
            .iwords = {bl.num_words, bl.detector_id, bl.subdet_id},
            .fwords = {abst0_val_},
            .td_refs = {tdr_id},
        });
        decode_partition_blocklet(dname, bl.data);
      }
    }
  }
  if (lraw_root <= 0) {
    if (auto it = bank_heads.find("RAWD"); it != bank_heads.end()) {
      for (const auto& b : ev.banks) {
        if (b.idh != "RAWD" || b.orig_addr <= 0) {
          continue;
        }
        // If a RAWD bank was not reached via RAW -> PMB, recover its partition from RAWD payload/ID
        process_rawd_chain("RAWD", 0, it->second);
        break;
      }
    }
  }

  // 3. Extract PV beamspot & PA particle candidates (pxdst33.car layout)
  reco::tracking::BeamspotConstraint beamspot{};
  struct ChargedCandidate {
    dst::TrackPidInput pid_input;
    reco::tracking::CandidateTrack solver_cand;
  };
  std::vector<ChargedCandidate> charged_candidates;
  std::int32_t next_syn_id = -1;

  auto ingest_pa = [&](std::size_t nd,
                       std::int32_t nmod,
                       std::int32_t iword0,
                       auto&& get_fword,
                       bool from_dtop_tree) {
    if (nd < 20 || nmod < 1 || nmod > 32) {
      return;
    }
    const std::size_t ip1 = static_cast<std::size_t>(4 + nmod);
    if (ip1 + 10 >= nd) {
      return;
    }
    const float px = get_fword(ip1 + 2);
    const float py = get_fword(ip1 + 3);
    const float pz = get_fword(ip1 + 4);
    const float energy = get_fword(ip1 + 5);
    const float ch_code_f = get_fword(ip1 + 7);
    const float rp = (ip1 + 10 < nd) ? get_fword(ip1 + 10) : 0.0f;
    const float dz = (ip1 + 19 < nd) ? get_fword(ip1 + 19) : 0.0f;

    if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(energy) || energy < 0.0f || energy > 200.0f) {
      return;
    }

    const int ch_code = static_cast<int>(std::lround(ch_code_f));
    if (ch_code == 1 || ch_code == 2) {
      const float pt = std::hypot(px, py);
      const float pmag = std::hypot(pt, pz);
      if (pt < 1e-3f || pmag > 500.0f) {
        return;
      }
      const std::int32_t q = (ch_code == 1) ? +1 : -1;
      const auto trk = track_from_p3_charge(
          px, py, pz, q, b0.bz,
          std::isfinite(rp) ? rp : 0.0f,
          std::isfinite(dz) ? dz : 0.0f);
      const std::int32_t tk_id =
          (!from_dtop_tree && iword0 > 0) ? iword0 : next_syn_id--;
      const double cand_score = static_cast<double>(energy) + 1.0;

      charged_candidates.push_back(ChargedCandidate{
          .pid_input =
              dst::TrackPidInput{
                  .track = trk,
                  .charge = q,
              },
          .solver_cand =
              reco::tracking::CandidateTrack{
                  .id = static_cast<std::int32_t>(charged_candidates.size()),
                  .te_ids = {tk_id},
                  .score = cand_score,
              },
      });
    } else {
      const float p_mag = std::sqrt(px * px + py * py + pz * pz);
      if (p_mag > 1e-4f && energy > 0.0f) {
        const float scale = 200.0f / p_mag;
        neutrals.push_back(dst::NeutralInput{
            .xyz_cm = {px * scale, py * scale, pz * scale},
            .energy_gev = energy,
            .is_ecal = true,
        });
      }
    }
  };

  bool used_dtop_tree = false;
  bool has_beamspot = false;
  if (auto it = bank_heads.find("DTOP"); it != bank_heads.end()) {
    const std::int32_t ldtop = it->second;
    if (zview.bank_nl(ldtop) >= 2 && zview.bank_link(ldtop, 2) > 0) {
      used_dtop_tree = true;
      std::int32_t lpv = zview.bank_link(ldtop, 2);
      std::size_t guard_pv = 0;
      while (lpv > 0 && guard_pv++ < 10'000) {
        // pxdst33.car:3277: skip simulated PV banks (JBIT(IQ(LPV), 4) != 0)
        if ((zview.bank_status(lpv) & 8) == 0) {
          if (!has_beamspot && zview.bank_nd(lpv) >= 7) {
            const float vx = zview.bank_fdata(lpv, 5);
            const float vy = zview.bank_fdata(lpv, 6);
            const float vz = zview.bank_fdata(lpv, 7);
            if (std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) &&
                std::fabs(vx) < 10.0f && std::fabs(vy) < 10.0f &&
                std::fabs(vz) < 50.0f) {
              beamspot.xyz_cm = {vx, vy, vz};
              has_beamspot = true;
            }
          }
          if (zview.bank_nl(lpv) >= 1) {
            std::int32_t lpa = zview.bank_link(lpv, 1);
            std::size_t guard_pa = 0;
            while (lpa > 0 && guard_pa++ < 10'000) {
              const std::int32_t nd_i = zview.bank_nd(lpa);
              if (nd_i >= 20) {
                const std::size_t nd = static_cast<std::size_t>(nd_i);
                const std::int32_t iword0 = zview.bank_idata(lpa, 1);
                const std::int32_t nmod = zview.bank_idata(lpa, 4);
                ingest_pa(
                    nd, nmod, iword0,
                    [&](std::size_t idx) {
                      return zview.bank_fdata(lpa, static_cast<std::int32_t>(idx) + 1);
                    },
                    true);
              }
              lpa = zview.bank_next(lpa);
            }
          }
        }
        lpv = zview.bank_next(lpv);
      }
    }
  }

  std::vector<const io::ParsedZebraBank*> pa_banks;
  if (!used_dtop_tree) {
    if (const auto* pv = ev.find_first("PV  ");
        pv != nullptr && pv->fwords.size() >= 7) {
      const float vx = pv->fwords[4];
      const float vy = pv->fwords[5];
      const float vz = pv->fwords[6];
      if (std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) &&
          std::fabs(vx) < 10.0f && std::fabs(vy) < 10.0f &&
          std::fabs(vz) < 50.0f) {
        beamspot.xyz_cm = {vx, vy, vz};
      }
    } else if (const auto* sp = ev.find_first("SP  ");
               sp != nullptr && sp->fwords.size() >= 4) {
      const float vx = sp->fwords[1];
      const float vy = sp->fwords[2];
      const float vz = sp->fwords[3];
      if (std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) &&
          std::fabs(vx) < 10.0f && std::fabs(vy) < 10.0f &&
          std::fabs(vz) < 50.0f) {
        beamspot.xyz_cm = {vx, vy, vz};
      }
    }

    pa_banks = ev.find_all("PA  ");
    for (const auto* pa : pa_banks) {
      ingest_pa(
          pa->fwords.size(),
          (pa->iwords.size() >= 4) ? pa->iwords[3] : 0,
          pa->iwords.empty() ? 0 : pa->iwords[0],
          [&](std::size_t idx) { return pa->fwords[idx]; },
          false);
    }
  }

  // Fallback for DELSIM MC truth SIM_ records: decode ST truth banks when no PA banks exist
  const auto st_banks = ev.find_all("ST  ");
  if (!used_dtop_tree && pa_banks.empty() && !st_banks.empty()) {
    for (std::size_t i = 0; i < st_banks.size(); ++i) {
      const auto* st = st_banks[i];
      if (st->fwords.size() < 8) {
        continue;
      }
      const float px = st->fwords[1];
      const float py = st->fwords[2];
      const float pz = st->fwords[3];
      const float energy = st->fwords[4];
      const float charge_f = st->fwords[6];
      const float mass = st->fwords[7];
      if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
          !std::isfinite(energy) || energy <= 0.01f || energy > 200.0f) {
        continue;
      }
      const int q = static_cast<int>(std::lround(charge_f));
      if (q == 1 || q == -1) {
        if (std::hypot(px, py) < 1e-3f) {
          continue;
        }
        const auto trk = track_from_p3_charge(px, py, pz, q, b0.bz);
        charged_candidates.push_back(ChargedCandidate{
            .pid_input = dst::TrackPidInput{.track = trk, .charge = q},
            .solver_cand = reco::tracking::CandidateTrack{
                .id = static_cast<std::int32_t>(charged_candidates.size()),
                .te_ids = {next_syn_id--},
                .score = static_cast<double>(energy) + 1.0,
            },
        });
      } else {
        const float p_mag = std::sqrt(px * px + py * py + pz * pz);
        if (p_mag > 1e-4f) {
          const float scale = 200.0f / p_mag;
          neutrals.push_back(dst::NeutralInput{
              .xyz_cm = {px * scale, py * scale, pz * scale},
              .energy_gev = energy,
              .is_ecal = (mass < 0.05f),
          });
        }
      }
    }
  }

  // 3b. Standalone pure-raw pattern recognition, Kalman refit & PID when no PA/ST exist
  if (!used_dtop_tree && pa_banks.empty() && st_banks.empty()) {
    std::sort(raw_track_hits.begin(), raw_track_hits.end(),
              [](const RawTrackHit& a, const RawTrackHit& b) {
                return a.meas.r_cm < b.meas.r_cm;
              });

    struct TrackCluster {
      float mean_phi = 0.0f;
      std::vector<const RawTrackHit*> hits;
    };
    std::vector<TrackCluster> clusters;
    for (const auto& h : raw_track_hits) {
      bool placed = false;
      for (auto& cl : clusters) {
        if (cl.hits.size() < 16u &&
            std::fabs(wrap_dphi(h.meas.phi_rad, cl.mean_phi)) < 0.28f) {
          const float w = static_cast<float>(cl.hits.size());
          cl.mean_phi = norm_phi(
              cl.mean_phi + wrap_dphi(h.meas.phi_rad, cl.mean_phi) / (w + 1.0f));
          cl.hits.push_back(&h);
          placed = true;
          break;
        }
      }
      if (!placed && clusters.size() < 32u) {
        clusters.push_back(TrackCluster{.mean_phi = h.meas.phi_rad, .hits = {&h}});
      }
    }

    for (const auto& cl : clusters) {
      if (cl.hits.size() < 2u) {
        continue;
      }
      const auto* h_in = cl.hits.front();
      const auto* h_out = cl.hits.back();
      const float dr = h_out->meas.r_cm - h_in->meas.r_cm;
      if (dr < 5.0f) {
        continue;
      }
      const float dphi_bend = wrap_dphi(h_out->meas.phi_rad, h_in->meas.phi_rad);
      float inv_r_seed = std::clamp(2.0f * dphi_bend / dr, -0.008f, 0.008f);
      if (std::fabs(inv_r_seed) < 1e-4f) {
        inv_r_seed = 0.0005f;
      }
      const float phi_seed = norm_phi(
          h_in->meas.phi_rad - 0.5f * inv_r_seed * h_in->meas.r_cm);

      double sum_th = 0.0;
      std::size_t n_th = 0;
      std::vector<reco::tracking::CylinderMeasurement> meas;
      std::vector<std::int32_t> te_ids;
      std::vector<float> dedx_samples;
      meas.reserve(cl.hits.size());
      te_ids.reserve(cl.hits.size());
      for (const auto* hp : cl.hits) {
        meas.push_back(hp->meas);
        te_ids.push_back(hp->ter_id > 0 ? hp->ter_id : next_syn_id--);
        if (hp->dedx_sample > 0.0f) {
          dedx_samples.push_back(hp->dedx_sample);
        }
        if (hp->meas.has_z && hp->meas.r_cm > 1.0f) {
          sum_th += std::atan2(static_cast<double>(hp->meas.r_cm),
                               static_cast<double>(hp->meas.z_cm));
          ++n_th;
        }
      }
      const float th_seed = (n_th > 0)
                                ? std::clamp(static_cast<float>(sum_th / static_cast<double>(n_th)),
                                             0.15f, kPi - 0.15f)
                                : (0.5f * kPi);

      const PerigeeTrack seed{
          .eps = 0.0f,
          .z0 = 0.0f,
          .theta = th_seed,
          .phi0 = phi_seed,
          .inv_r = inv_r_seed,
      };
      const auto fit = reco::tracking::kalman_refit(seed, b0.bz, meas);
      ++out.kalman_refits;

      PerigeeTrack chosen = fit.track;
      const float pt_fit = chosen.pt_gev(b0.bz);
      const float p_fit = chosen.p_gev(b0.bz);
      if (!std::isfinite(chosen.eps) || !std::isfinite(chosen.z0) ||
          !std::isfinite(chosen.theta) || !std::isfinite(chosen.phi0) ||
          !std::isfinite(chosen.inv_r) || !std::isfinite(pt_fit) ||
          !std::isfinite(p_fit) || pt_fit < 0.2f || pt_fit > 60.0f ||
          p_fit > 80.0f || chosen.theta < 0.05f || chosen.theta > kPi - 0.05f) {
        chosen = seed;
      }

      const std::int32_t q = (chosen.inv_r >= 0.0f) ? +1 : -1;
      const float trk_phi = norm_phi(chosen.phi0);
      float ecal_e = 0.0f;
      float hcal_e = 0.0f;
      for (auto& cc : raw_calo_clusters) {
        if (std::fabs(wrap_dphi(cc.phi_rad, trk_phi)) < 0.25f &&
            std::fabs(cc.theta_rad - chosen.theta) < 0.35f) {
          if (cc.is_ecal) {
            ecal_e += cc.energy_gev;
          } else {
            hcal_e += cc.energy_gev;
          }
          cc.matched = true;
        }
      }

      bool has_mu = false;
      for (const auto& mh : raw_muon_hits) {
        if (std::fabs(wrap_dphi(mh.phi_rad, trk_phi)) < 0.30f) {
          has_mu = true;
          break;
        }
      }

      bool has_rich = false;
      reco::rich::Species rich_sp = reco::rich::Species::Pion;
      {
        double sum_tc = 0.0;
        std::int32_t n_ph = 0;
        std::int32_t rich_det = reco::rich::kDetectorIdBrich;
        for (const auto& rh : raw_rich_hits) {
          if (std::fabs(wrap_dphi(rh.phi_rad, trk_phi)) < 0.25f) {
            sum_tc += static_cast<double>(rh.theta_c_rad);
            ++n_ph;
            rich_det = rh.detector_id;
          }
        }
        if (n_ph > 0) {
          const float p_mag = std::min(chosen.p_gev(b0.bz), 45.0f);
          const reco::rich::RichRingMeasurement ring{
              .momentum_gev = p_mag,
              .refractive_index = reco::rich::kRefractiveIndexLiquid,
              .measured_theta_rad = static_cast<float>(sum_tc / static_cast<double>(n_ph)),
              .sigma_theta_rad = 0.01f,
              .num_photons = n_ph,
          };
          const auto lr = reco::rich::evaluate_hypotheses(ring);
          reco::rich::write_ter_bank(out.tanagra, ring, lr, rich_det);
          has_rich = true;
          rich_sp = lr.best_species;
        }
      }

      const float dedx = !dedx_samples.empty()
                             ? reco::tpc::truncated_mean_dedx(dedx_samples)
                             : 1.0f;

      charged_candidates.push_back(ChargedCandidate{
          .pid_input = dst::TrackPidInput{
              .track = chosen,
              .charge = q,
              .ecal_energy_gev = ecal_e,
              .hcal_energy_gev = hcal_e,
              .tpc_dedx = dedx,
              .has_muon_hit = has_mu,
              .rich_species = rich_sp,
              .has_rich = has_rich,
          },
          .solver_cand = reco::tracking::CandidateTrack{
              .id = static_cast<std::int32_t>(charged_candidates.size()),
              .te_ids = std::move(te_ids),
              .score = static_cast<double>(cl.hits.size()) +
                       static_cast<double>(ecal_e + hcal_e) * 0.1,
          },
      });
    }

    for (const auto& cc : raw_calo_clusters) {
      if (!cc.matched) {
        neutrals.push_back(dst::NeutralInput{
            .xyz_cm = cc.xyz_cm,
            .energy_gev = cc.energy_gev,
            .is_ecal = cc.is_ecal,
        });
      }
    }
  }

  // 4. Global ambiguity resolution & TKR banking
  out.tracks_before_ambiguity = charged_candidates.size();
  std::vector<dst::TrackPidInput> accepted_tracks;

  if (!charged_candidates.empty()) {
    std::vector<reco::tracking::CandidateTrack> solver_in;
    solver_in.reserve(charged_candidates.size());
    for (const auto& c : charged_candidates) {
      solver_in.push_back(c.solver_cand);
    }
    const auto solved = reco::tracking::solve_ambiguities(solver_in);
    std::unordered_set<std::int32_t> keep(
        solved.accepted_track_ids.begin(), solved.accepted_track_ids.end());

    for (const auto& c : charged_candidates) {
      if (!keep.contains(c.solver_cand.id)) {
        continue;
      }
      accepted_tracks.push_back(c.pid_input);
      const auto& trk = c.pid_input.track;
      out.tanagra.add_tkr(TkrBank{
          .charge = c.pid_input.charge,
          .perigee_params = {trk.eps, trk.z0, trk.theta, trk.phi0, trk.inv_r},
      });
    }
  }
  out.tracks_after_ambiguity = accepted_tracks.size();

  // 5. Short-DST construction (primary vertex + PID + kinematics)
  out.dst_event = dst::build_short_dst(
      out.run, out.event, b0.bz, accepted_tracks, neutrals, beamspot, is_bhabha);

  return out;
}

}  // namespace delphipp::pipeline
