// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/pipeline/event_pipeline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_set>
#include <vector>

#include "delphipp/helix.hpp"
#include "delphipp/io/raw_blocklet.hpp"
#include "delphipp/reco/ambiguity_solver.hpp"
#include "delphipp/reco/hof.hpp"
#include "delphipp/reco/kalman_fitter.hpp"
#include "delphipp/reco/vertex_fitter.hpp"
#include "delphipp/zebra_store.hpp"

namespace delphipp::pipeline {

namespace {

PerigeeTrack track_from_p3_charge(float px,
                                  float py,
                                  float pz,
                                  std::int32_t charge,
                                  float bz,
                                  float eps = 0.0f,
                                  float z0 = 0.0f) {
  const float p = std::sqrt(px * px + py * py + pz * pz);
  const float theta =
      (p > 1e-6f) ? std::acos(std::clamp(pz / p, -1.0f, 1.0f)) : 1.5707963f;
  const float phi0 = std::atan2(py, px);
  const float signed_q = (charge == 0) ? 1.0f : static_cast<float>(charge);
  const float inv_r = (p > 1e-5f) ? (signed_q * bz / p) : 0.0f;
  return PerigeeTrack{
      .eps = eps,
      .z0 = z0,
      .theta = theta,
      .phi0 = phi0,
      .inv_r = inv_r,
  };
}

constexpr std::array<const char*, 18> kRawBankIds = {
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
  out.run = (ev.run_number > 0) ? ev.run_number : fallback_run;
  out.event = (ev.event_number > 0) ? ev.event_number : fallback_evt;
  out.zebra_banks = ev.banks.size();

  const auto b0 = bfield_.evaluate({0.0f, 0.0f, 0.0f});

  OwnedZebraStore zstore(65536);
  auto zview = zstore.view();
  const auto bank_heads = ev.populate_zebra_store(zview);

  // 1. Specialized HOF raw-to-TD unpacking when EGRA/EQRA banks are present
  reco::hof::HofCalibration hof_calib{};
  reco::hof::HofStats hof_stats{};
  hof_stats.begin_event();
  if (auto it = bank_heads.find("EGRA"); it != bank_heads.end()) {
    const std::int32_t mraw =
        bank_heads.contains("EQRA") ? bank_heads.at("EQRA") : 0;
    const auto hof_res =
        reco::hof::unpack_td(zview, it->second, mraw, hof_calib, hof_stats);
    out.hof_hits = hof_res.counter_hits;
    if (!hof_res.hits.empty()) {
      std::vector<float> flat;
      flat.reserve(hof_res.hits.size() * 8);
      for (const auto& h : hof_res.hits) {
        flat.insert(flat.end(), h.words.begin(), h.words.end());
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

  // 2. Uniform raw blocklet decoding across all 18 DELPHI detector banks
  for (const char* raw_id : kRawBankIds) {
    if (auto it = bank_heads.find(raw_id); it != bank_heads.end()) {
      const auto blocklets = io::read_raw_blocklet_chain(zview, it->second);
      out.raw_blocklets += blocklets.size();
      for (const auto& bl : blocklets) {
        std::string dname = detector_name(bl.detector_id);
        if (std::string(raw_id) == "EGRA" && bl.detector_id == reco::hof::kModuleSim) {
          dname = "HOF";
        }
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
      }
    }
  }

  // 3. Extract PV beamspot & PA particle candidates (pxdst33.car layout)
  reco::tracking::BeamspotConstraint beamspot{};
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
  }

  struct ChargedCandidate {
    dst::TrackPidInput pid_input;
    reco::tracking::CandidateTrack solver_cand;
  };
  std::vector<ChargedCandidate> charged_candidates;
  std::vector<dst::NeutralInput> neutrals;

  std::int32_t next_syn_id = -1;
  const auto pa_banks = ev.find_all("PA  ");
  for (std::size_t i = 0; i < pa_banks.size(); ++i) {
    const auto* pa = pa_banks[i];
    if (pa->iwords.size() < 20) {
      continue;
    }
    const std::int32_t nmod = pa->iwords[3];
    if (nmod < 1 || nmod > 32) {
      continue;
    }
    const std::size_t ip1 = static_cast<std::size_t>(4 + nmod);
    if (ip1 + 10 >= pa->fwords.size()) {
      continue;
    }
    const float px = pa->fwords[ip1 + 2];
    const float py = pa->fwords[ip1 + 3];
    const float pz = pa->fwords[ip1 + 4];
    const float energy = pa->fwords[ip1 + 5];
    const float ch_code_f = pa->fwords[ip1 + 7];
    const float rp =
        (ip1 + 10 < pa->fwords.size()) ? pa->fwords[ip1 + 10] : 0.0f;
    const float dz =
        (ip1 + 19 < pa->fwords.size()) ? pa->fwords[ip1 + 19] : 0.0f;

    if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(energy) || energy < 0.0f || energy > 200.0f) {
      continue;
    }

    const int ch_code = static_cast<int>(std::lround(ch_code_f));
    if (ch_code == 1 || ch_code == 2) {
      const std::int32_t q = (ch_code == 1) ? +1 : -1;
      const auto trk = track_from_p3_charge(
          px, py, pz, q, b0.bz,
          std::isfinite(rp) ? rp : 0.0f,
          std::isfinite(dz) ? dz : 0.0f);
      const std::int32_t tk_id =
          (pa->iwords[0] > 0) ? pa->iwords[0] : next_syn_id--;
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
      if (p_mag > 1e-4f) {
        const float scale = 200.0f / p_mag;
        neutrals.push_back(dst::NeutralInput{
            .xyz_cm = {px * scale, py * scale, pz * scale},
            .energy_gev = energy,
            .is_ecal = true,
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
      out.run, out.event, b0.bz, accepted_tracks, neutrals, beamspot);

  return out;
}

}  // namespace delphipp::pipeline
