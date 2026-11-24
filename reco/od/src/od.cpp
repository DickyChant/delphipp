// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/od.hpp"

namespace delphipp::reco::od {

std::vector<OdTrackElement> reconstruct_planks(std::span<const OdTubeHit> hits,
                                               float v_prop_cm_per_ns) {
  std::vector<OdTrackElement> out;
  for (std::int32_t p = 1; p <= kNumPlanks; ++p) {
    double sum_r = 0.0, sum_rp = 0.0, sum_z = 0.0;
    std::int32_t n = 0;
    std::int32_t nz = 0;
    for (const auto& h : hits) {
      if (h.plank != p) continue;
      sum_r += static_cast<double>(h.r_cm);
      sum_rp += static_cast<double>(h.r_phi_cm);
      ++n;
      if (h.has_z_timing) {
        sum_z += 0.5 * static_cast<double>(v_prop_cm_per_ns) *
                 static_cast<double>(h.t_end_c_ns - h.t_end_a_ns);
        ++nz;
      }
    }
    if (n >= 2) {
      out.push_back(OdTrackElement{
          .plank = p,
          .r_cm = static_cast<float>(sum_r / static_cast<double>(n)),
          .r_phi_cm = static_cast<float>(sum_rp / static_cast<double>(n)),
          .z_cm = (nz > 0) ? static_cast<float>(sum_z / static_cast<double>(nz)) : 0.0f,
          .num_layers = n,
      });
    }
  }
  return out;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const OdTrackElement> elements) {
  std::vector<TanagraId> ids;
  ids.reserve(elements.size());
  for (const auto& e : elements) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = e.plank,
        .descriptor = 125,
        .fwords = {e.r_cm, e.r_phi_cm, e.z_cm},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::od
