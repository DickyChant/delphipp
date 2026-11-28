// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/trigger/trigger.hpp"

#include <algorithm>
#include <cstddef>

namespace delphipp::reco::trigger {

namespace {

// Exact T2 collinearity sum window from fcabtrig.car line 948:
//   PARAMETER (MNT2=93, MXT2=100)
inline constexpr std::int32_t kFcabT2MinSum = 93;
inline constexpr std::int32_t kFcabT2MaxSum = 100;

struct EndcapEval {
  std::int32_t t1 = 0;
  std::int32_t t2 = 0;
  std::int32_t n_3coord = 0;
};

[[nodiscard]] EndcapEval evaluate_fcab_endcap(const FcabEndcapInput& ec) noexcept {
  EndcapEval out{};
  // Majority logic in each FCB projection (V, U, X): >= 3 of 4 planes
  std::int32_t proj_trig = 0;
  for (std::size_t d = 0; d < 3; ++d) {
    if (ec.fcb_plane_hits[d] >= 3) {
      ++proj_trig;
    }
  }
  out.n_3coord = proj_trig;

  // T1 requires FCA presence and >= 2 of 3 FCB coordinate projections
  if (ec.fca_hit_groups > 0 && proj_trig >= 2) {
    out.t1 = std::min<std::int32_t>(proj_trig, 3);
  } else if (proj_trig == 3) {
    out.t1 = 1;
  }

  // T2S sum condition (fcabtrig.car lines 1111-1121)
  if (!ec.t1b_groups[0].empty() && !ec.t1b_groups[1].empty() &&
      !ec.t1b_groups[2].empty()) {
    std::int32_t n_cand = 0;
    const std::size_t nv = std::min<std::size_t>(16u, ec.t1b_groups[0].size());
    const std::size_t nu = std::min<std::size_t>(16u, ec.t1b_groups[1].size());
    const std::size_t nx = std::min<std::size_t>(16u, ec.t1b_groups[2].size());
    for (std::size_t iv = 0; iv < nv; ++iv) {
      for (std::size_t iu = 0; iu < nu; ++iu) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
          const std::int32_t sum =
              ec.t1b_groups[0][iv] + ec.t1b_groups[1][iu] + ec.t1b_groups[2][ix];
          if (sum >= kFcabT2MinSum && sum <= kFcabT2MaxSum) {
            ++n_cand;
            if (n_cand >= 3) {
              out.t2 = 3;
              return out;
            }
          }
        }
      }
    }
    out.t2 = n_cand;
  }

  return out;
}

}  // namespace

FcabTriggerResult evaluate_fcab_trigger(const FcabTriggerInput& input) noexcept {
  const auto fw = evaluate_fcab_endcap(input.forward);
  const auto bw = evaluate_fcab_endcap(input.backward);
  return FcabTriggerResult{
      .t1_forward = fw.t1,
      .t1_backward = bw.t1,
      .t2_forward = fw.t2,
      .t2_backward = bw.t2,
      .n_3coord_forward = fw.n_3coord,
      .n_3coord_backward = bw.n_3coord,
  };
}

TriggerDecision evaluate_trigger(const TriggerDetectorInput& input,
                                 const T1Mask& mask) noexcept {
  TriggerDecision dec{};
  dec.fcab = evaluate_fcab_trigger(input.fcab);

  const bool has_fcab_t1_fw = input.fca_fcb_fw || (dec.fcab.t1_forward > 0);
  const bool has_fcab_t1_bw = input.fca_fcb_bw || (dec.fcab.t1_backward > 0);
  const bool has_fcab_t2_fw = input.fcab_t2_fw || (dec.fcab.t2_forward > 0);
  const bool has_fcab_t2_bw = input.fcab_t2_bw || (dec.fcab.t2_backward > 0);

  // T1 subtrigger coincidences (trigger.car lines 93-97, 320-336)
  dec.t1.id_od = (input.id_hits >= 1 && input.od_hits >= 2);
  dec.t1.tof_od = (input.tof_hits >= 1 && input.od_hits >= 1);
  dec.t1.tof_bb = (input.tof_hits >= 2 && input.tof_back_to_back);
  dec.t1.hpc_od = (input.hpc_clusters >= 1 && input.od_hits >= 1);
  dec.t1.hpc_tof = (input.hpc_clusters >= 1 && input.tof_hits >= 1);
  dec.t1.hpc_bb = (input.hpc_clusters >= 2 && input.hpc_back_to_back);
  dec.t1.hpc_maj = (input.hpc_clusters >= 2) || (input.hpc_high_clusters >= 1);

  const bool fwd_a = has_fcab_t1_fw && (input.emf_low_fw || input.hof_hits_fw > 0);
  const bool fwd_c = has_fcab_t1_bw && (input.emf_low_bw || input.hof_hits_bw > 0);
  const bool fca_fcb_coinc = has_fcab_t1_fw && has_fcab_t1_bw;
  dec.t1.fwd_maj = fwd_a || fwd_c || fca_fcb_coinc;

  dec.t1.emf_bb =
      (input.emf_high_fw && input.emf_high_bw) ||
      (input.emf_low_fw && input.emf_low_bw) ||
      (input.emf_energy_fw_gev >= 3.0f && input.emf_energy_bw_gev >= 3.0f);

  dec.t1.mub_coinc =
      (input.mub_hits >= 1 &&
       (input.id_hits >= 1 || input.od_hits >= 1 || input.tof_hits >= 1)) ||
      (input.mub_hits >= 2 && input.mub_back_to_back);

  dec.t1.sat_bhabha = input.sat_bhabha_coincidence;

  dec.t1.t1_accept =
      (mask.id_od && dec.t1.id_od) ||
      (mask.tof_od && dec.t1.tof_od) ||
      (mask.tof_bb && dec.t1.tof_bb) ||
      (mask.hpc_od && dec.t1.hpc_od) ||
      (mask.hpc_tof && dec.t1.hpc_tof) ||
      (mask.hpc_bb && dec.t1.hpc_bb) ||
      (mask.hpc_maj && dec.t1.hpc_maj) ||
      (mask.fwd_maj && dec.t1.fwd_maj) ||
      (mask.emf_bb && dec.t1.emf_bb) ||
      (mask.mub_coinc && dec.t1.mub_coinc) ||
      (mask.sat_bhabha && dec.t1.sat_bhabha);

  if (!dec.t1.t1_accept) {
    return dec;
  }

  // T2 confirmation signals (trigger.car lines 155-159, 1466-1478)
  dec.t2.tpc_confirm = (input.tpc_tracks >= 2) || input.tpc_back_to_back;
  dec.t2.hpc_confirm =
      (input.hpc_clusters >= 1) &&
      (input.hpc_energy_gev >= 1.5f || input.hpc_high_clusters >= 1);
  dec.t2.hcal_confirm =
      (input.hcal_clusters >= 1 && input.hcal_energy_gev >= 2.0f);
  dec.t2.fcab_confirm = has_fcab_t2_fw || has_fcab_t2_bw;
  dec.t2.mub_confirm =
      (input.mub_hits >= 1) &&
      (input.tpc_tracks >= 1 || input.id_hits >= 1 || input.od_hits >= 1);
  dec.t2.emf_sat_confirm = dec.t1.emf_bb || dec.t1.sat_bhabha;

  dec.t2.t2_accept =
      dec.t2.tpc_confirm || dec.t2.hpc_confirm || dec.t2.hcal_confirm ||
      dec.t2.fcab_confirm || dec.t2.mub_confirm || dec.t2.emf_sat_confirm;

  if (!dec.t2.t2_accept) {
    return dec;
  }

  // T3 online event filter (trigger.car TrGetT3)
  const float total_em =
      input.hpc_energy_gev + input.emf_energy_fw_gev + input.emf_energy_bw_gev;
  const float total_calo = total_em + input.hcal_energy_gev;

  dec.t3_accept =
      (input.charged_multiplicity >= 1) ||
      (input.tpc_tracks >= 1) ||
      (total_calo >= 2.0f) ||
      dec.t1.sat_bhabha;

  if (!dec.t3_accept) {
    return dec;
  }

  // T4 physics stream classification: HadronicZ, Bhabha, Dilepton, GammaGamma
  const bool is_hadronic_z =
      (input.charged_multiplicity >= 5) &&
      (input.visible_energy_gev >= 15.0f || total_calo >= 12.0f);

  const bool is_barrel_bhabha =
      (input.charged_multiplicity >= 1 && input.charged_multiplicity <= 4) &&
      (input.hpc_energy_gev >= 30.0f) &&
      (input.hpc_back_to_back || dec.t1.hpc_bb || dec.t2.tpc_confirm) &&
      (input.hcal_energy_gev < 10.0f);
  const bool is_endcap_bhabha =
      (input.charged_multiplicity <= 4) &&
      (input.emf_energy_fw_gev >= 10.0f && input.emf_energy_bw_gev >= 10.0f) &&
      ((input.emf_energy_fw_gev + input.emf_energy_bw_gev) >= 30.0f) &&
      (input.hcal_energy_gev < 10.0f);
  const bool is_bhabha =
      is_barrel_bhabha || is_endcap_bhabha || dec.t1.sat_bhabha;

  const bool has_muon_sig = (input.mub_hits >= 1 || input.muf_hits >= 1);
  const bool is_dilepton =
      (!is_bhabha || has_muon_sig) &&
      (input.charged_multiplicity >= 1 && input.charged_multiplicity <= 4) &&
      (input.max_track_p_gev >= 10.0f || input.visible_energy_gev >= 25.0f) &&
      (has_muon_sig || dec.t1.tof_bb || dec.t2.tpc_confirm ||
       input.hpc_energy_gev >= 8.0f);

  const bool is_gamma_gamma =
      !is_hadronic_z && !is_bhabha && !is_dilepton &&
      (input.charged_multiplicity <= 6) &&
      ((input.visible_energy_gev >= 1.5f && input.visible_energy_gev < 20.0f) ||
       (total_em >= 1.0f && total_em < 15.0f &&
        (input.hpc_clusters >= 2 || (input.emf_low_fw && input.emf_low_bw) ||
         input.charged_multiplicity >= 2)));

  if (is_hadronic_z) {
    dec.streams = dec.streams | PhysicsStream::HadronicZ;
  }
  if (is_bhabha) {
    dec.streams = dec.streams | PhysicsStream::Bhabha;
  }
  if (is_dilepton) {
    dec.streams = dec.streams | PhysicsStream::Dilepton;
  }
  if (is_gamma_gamma) {
    dec.streams = dec.streams | PhysicsStream::GammaGamma;
  }

  // Primary stream priority: HadronicZ > Bhabha > Dilepton > GammaGamma
  if (is_hadronic_z) {
    dec.primary_stream = PhysicsStream::HadronicZ;
  } else if (is_bhabha) {
    dec.primary_stream = PhysicsStream::Bhabha;
  } else if (is_dilepton) {
    dec.primary_stream = PhysicsStream::Dilepton;
  } else if (is_gamma_gamma) {
    dec.primary_stream = PhysicsStream::GammaGamma;
  }

  dec.t4_accept = (dec.streams != PhysicsStream::None);
  return dec;
}

}  // namespace delphipp::reco::trigger
