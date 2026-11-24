// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/id.hpp"

#include <cmath>
#include <numbers>

namespace delphipp::reco::id {

std::vector<IdSegment> reconstruct_segments(std::span<const IdJetHit> hits,
                                            float v_drift_cm_per_ns) {
  constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
  constexpr float kSectorWidth = kTwoPi / static_cast<float>(kNumSectors);
  std::vector<IdSegment> segs;

  for (std::int32_t sec = 1; sec <= kNumSectors; ++sec) {
    const float sec_phi = (static_cast<float>(sec - 1) + 0.5f) * kSectorWidth;
    double sr = 0.0, sp = 0.0, srr = 0.0, srp = 0.0;
    std::int32_t n = 0;

    for (const auto& h : hits) {
      if (h.sector != sec || h.r_cm <= 0.0f) continue;
      const float drift_cm = static_cast<float>(h.lr_sign) * v_drift_cm_per_ns * h.drift_time_ns;
      const double phi = static_cast<double>(sec_phi + drift_cm / h.r_cm);
      const double r = static_cast<double>(h.r_cm);
      sr += r;
      sp += phi;
      srr += r * r;
      srp += r * phi;
      ++n;
    }

    if (n >= 2) {
      const double mean_r = sr / static_cast<double>(n);
      const double mean_p = sp / static_cast<double>(n);
      const double denom = srr - static_cast<double>(n) * mean_r * mean_r;
      const double slope = (std::fabs(denom) > 1e-9)
                               ? (srp - static_cast<double>(n) * mean_r * mean_p) / denom
                               : 0.0;
      segs.push_back(IdSegment{
          .sector = sec,
          .r_ref_cm = static_cast<float>(mean_r),
          .phi_rad = static_cast<float>(mean_p),
          .dphi_dr = static_cast<float>(slope),
          .num_hits = n,
      });
    }
  }
  return segs;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const IdSegment> segments) {
  std::vector<TanagraId> ids;
  ids.reserve(segments.size());
  for (const auto& s : segments) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = s.sector,
        .descriptor = 125,
        .fwords = {s.r_ref_cm, s.r_ref_cm * s.phi_rad, 0.0f, s.dphi_dr},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::id
