// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/fwd_chambers.hpp"

#include <array>
#include <cmath>
#include <numbers>

namespace delphipp::reco::fwd {

std::vector<ForwardWireHit> unpack_raw_blocklet(std::span<const std::int32_t> words,
                                                std::int32_t detector_id) {
  std::vector<ForwardWireHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  constexpr float kPi = std::numbers::pi_v<float>;
  constexpr std::array<float, 3> kAngles = {0.0f, kPi / 3.0f, -kPi / 3.0f};
  const float abs_z = (detector_id == kDetectorIdFcb) ? 275.0f : 160.0f;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::uint32_t u_raw = u & 0x3FFu;
    if (u_raw == 0u) {
      continue;
    }
    const std::size_t plane = static_cast<std::size_t>((u >> 10) & 0x3u) % 3u;
    const std::int32_t endcap = ((u >> 12) & 1u) ? 2 : 1;
    hits.push_back(ForwardWireHit{
        .endcap = endcap,
        .z_cm = (endcap == 1) ? abs_z : -abs_z,
        .alpha_rad = kAngles[plane],
        .u_cm = (static_cast<float>(u_raw) - 512.0f) * 0.25f,
        .sigma_u_cm = 0.03f,
    });
  }
  return hits;
}

std::optional<ForwardSegment> fit_segment(std::span<const ForwardWireHit> hits) noexcept {
  if (hits.size() < 2) {
    return std::nullopt;
  }
  double a11 = 0.0, a12 = 0.0, a22 = 0.0;
  double b1 = 0.0, b2 = 0.0;
  double sum_z = 0.0;

  for (const auto& h : hits) {
    const double w = 1.0 / (static_cast<double>(h.sigma_u_cm) * static_cast<double>(h.sigma_u_cm));
    const double c = std::cos(static_cast<double>(h.alpha_rad));
    const double s = std::sin(static_cast<double>(h.alpha_rad));
    const double u = static_cast<double>(h.u_cm);
    a11 += w * c * c;
    a12 += w * c * s;
    a22 += w * s * s;
    b1 += w * c * u;
    b2 += w * s * u;
    sum_z += static_cast<double>(h.z_cm);
  }

  const double det = a11 * a22 - a12 * a12;
  if (std::fabs(det) < 1e-9) {
    return std::nullopt;
  }

  const double x = (a22 * b1 - a12 * b2) / det;
  const double y = (-a12 * b1 + a11 * b2) / det;

  double chi2 = 0.0;
  for (const auto& h : hits) {
    const double c = std::cos(static_cast<double>(h.alpha_rad));
    const double s = std::sin(static_cast<double>(h.alpha_rad));
    const double res = (x * c + y * s - static_cast<double>(h.u_cm)) /
                       static_cast<double>(h.sigma_u_cm);
    chi2 += res * res;
  }

  return ForwardSegment{
      .endcap = hits[0].endcap,
      .x_cm = static_cast<float>(x),
      .y_cm = static_cast<float>(y),
      .z_cm = static_cast<float>(sum_z / static_cast<double>(hits.size())),
      .chi2 = static_cast<float>(chi2),
      .num_planes = static_cast<std::int32_t>(hits.size()),
  };
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const ForwardSegment> segments,
                                       std::int32_t detector_id) {
  std::vector<TanagraId> ids;
  ids.reserve(segments.size());
  for (const auto& seg : segments) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = detector_id,
        .submod_id = seg.endcap,
        .descriptor = 124,
        .fwords = {seg.x_cm, seg.y_cm, seg.z_cm, seg.chi2},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::fwd
