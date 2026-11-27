// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/stw.hpp"

#include <cmath>
#include <numbers>

namespace delphipp::reco::stw {

std::vector<StwHit> unpack_raw_blocklet(std::span<const std::int32_t> words) {
  std::vector<StwHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;

  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::uint32_t tdc = u & 0xFFu;
    if (tdc == 0u) {
      continue;
    }
    const std::int32_t straw = static_cast<std::int32_t>((u >> 8) & 0xFFu) % kNumStrawsPerLayer + 1;
    const std::int32_t layer = static_cast<std::int32_t>((u >> 16) & 0x7u) % kNumLayers + 1;
    const std::int32_t sector = (straw - 1) / (kNumStrawsPerLayer / kNumSectors) + 1;
    const float r_cm = kLayerRadiiCm[static_cast<std::size_t>(layer - 1)];
    const float phi_rad = (static_cast<float>(straw - 1) + 0.5f) *
                          (kTwoPi / static_cast<float>(kNumStrawsPerLayer));
    hits.push_back(StwHit{
        .layer = layer,
        .straw = straw,
        .sector = sector,
        .r_cm = r_cm,
        .phi_rad = phi_rad,
        .drift_time_ns = static_cast<float>(tdc) * 2.0f,
        .pulse_height = static_cast<float>((u >> 19) & 0xFFu),
    });
  }
  return hits;
}

std::vector<StwTrackElement> reconstruct_elements(std::span<const StwHit> hits,
                                                  float v_drift_cm_per_ns) {
  constexpr float kPi = std::numbers::pi_v<float>;
  constexpr float kTwoPi = 2.0f * kPi;
  std::vector<StwTrackElement> elems;

  for (std::int32_t sec = 1; sec <= kNumSectors; ++sec) {
    double sr = 0.0, sp = 0.0, srr = 0.0, srp = 0.0;
    double phi_ref = 0.0;
    bool has_ref = false;
    std::int32_t n = 0;

    for (const auto& h : hits) {
      if (h.sector != sec || h.r_cm <= 0.0f) {
        continue;
      }
      const float drift_cm = v_drift_cm_per_ns * h.drift_time_ns;
      double phi = static_cast<double>(h.phi_rad) +
                   static_cast<double>(drift_cm / h.r_cm);
      if (!has_ref) {
        phi_ref = phi;
        has_ref = true;
      } else {
        while (phi - phi_ref > static_cast<double>(kPi)) phi -= static_cast<double>(kTwoPi);
        while (phi - phi_ref < -static_cast<double>(kPi)) phi += static_cast<double>(kTwoPi);
      }
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
      double norm_p = mean_p;
      while (norm_p < 0.0) norm_p += static_cast<double>(kTwoPi);
      while (norm_p >= static_cast<double>(kTwoPi)) norm_p -= static_cast<double>(kTwoPi);

      elems.push_back(StwTrackElement{
          .sector = sec,
          .r_ref_cm = static_cast<float>(mean_r),
          .phi_rad = static_cast<float>(norm_p),
          .dphi_dr = static_cast<float>(slope),
          .num_layers = n,
      });
    }
  }
  return elems;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const StwTrackElement> elements) {
  std::vector<TanagraId> ids;
  ids.reserve(elements.size());
  for (const auto& e : elements) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = e.sector,
        .descriptor = 125,
        .fwords = {e.r_ref_cm, e.r_ref_cm * e.phi_rad, 0.0f, e.dphi_dr},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::stw
