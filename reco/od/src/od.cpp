// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/od.hpp"

#include <numbers>

namespace delphipp::reco::od {

std::vector<OdTubeHit> unpack_raw_blocklet(std::span<const std::int32_t> words) {
  std::vector<OdTubeHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::uint32_t tdc_a = u & 0xFFu;
    const std::uint32_t tdc_c = (u >> 8) & 0xFFu;
    if (tdc_a == 0u && tdc_c == 0u) {
      continue;
    }
    const std::int32_t layer = static_cast<std::int32_t>((u >> 16) & 0x7u) % 5 + 1;
    const std::int32_t plank = static_cast<std::int32_t>((u >> 19) & 0x1Fu) % kNumPlanks + 1;
    const float r_cm = 197.0f + static_cast<float>(layer - 1) * 2.25f;
    const float phi = (static_cast<float>(plank - 1) + 0.5f) * (kTwoPi / static_cast<float>(kNumPlanks));
    hits.push_back(OdTubeHit{
        .plank = plank,
        .layer = layer,
        .r_cm = r_cm,
        .r_phi_cm = r_cm * phi,
        .t_end_a_ns = static_cast<float>(tdc_a) * 0.5f,
        .t_end_c_ns = static_cast<float>(tdc_c) * 0.5f,
        .has_z_timing = (tdc_a > 0u && tdc_c > 0u),
    });
  }
  return hits;
}

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
