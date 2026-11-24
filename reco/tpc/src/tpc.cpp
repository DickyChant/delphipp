// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/tpc.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace delphipp::reco::tpc {

float reconstruct_z(std::int32_t endplate,
                    float time_bin,
                    const TpcCalibration& calib) noexcept {
  const float drift_us = (time_bin - calib.abst0_bins) * calib.time_bin_us;
  const float drift_dist = std::max(0.0f, calib.drift_velocity_cm_per_us * drift_us);
  const float abs_z = calib.endplate_z_cm - drift_dist;
  return (endplate == 0) ? abs_z : -abs_z;
}

std::vector<TpcSpacePoint> reconstruct_space_points(
    std::span<const TpcPadHit> pad_hits,
    const TpcCalibration& calib) {
  std::vector<TpcSpacePoint> pts;
  pts.reserve(pad_hits.size());
  for (const auto& h : pad_hits) {
    pts.push_back(TpcSpacePoint{
        .endplate = h.endplate,
        .sector = h.sector,
        .pad_row = h.pad_row,
        .r_cm = h.r_cm,
        .r_phi_cm = h.r_phi_cm,
        .z_cm = reconstruct_z(h.endplate, h.time_bin, calib),
        .charge = h.charge,
    });
  }
  return pts;
}

float truncated_mean_dedx(std::span<const float> wire_samples,
                          float keep_fraction) {
  if (wire_samples.empty()) {
    return 0.0f;
  }
  std::vector<float> sorted;
  sorted.reserve(wire_samples.size());
  for (const float s : wire_samples) {
    if (s > 0.0f) {
      sorted.push_back(s);
    }
  }
  if (sorted.empty()) {
    return 0.0f;
  }
  std::sort(sorted.begin(), sorted.end());
  const std::size_t n_keep = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor(static_cast<float>(sorted.size()) * keep_fraction)));
  double sum = 0.0;
  for (std::size_t i = 0; i < n_keep; ++i) {
    sum += static_cast<double>(sorted[i]);
  }
  return static_cast<float>(sum / static_cast<double>(n_keep));
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const TpcSpacePoint> points) {
  std::vector<TanagraId> ids;
  ids.reserve(points.size());
  for (const auto& p : points) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = p.endplate * 100 + p.sector,
        .descriptor = 125,
        .fwords = {p.r_cm, p.r_phi_cm, p.z_cm, p.charge},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::tpc
