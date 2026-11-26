// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/sat.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace delphipp::reco::sat {

std::vector<SatElementHit> unpack_raw_blocklet(std::span<const std::int32_t> words) {
  std::vector<SatElementHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::uint32_t adc = u & 0xFFFu;
    if (adc == 0u) {
      continue;
    }
    const std::int32_t ring = static_cast<std::int32_t>((u >> 12) & 0x7u) + 1;
    const std::int32_t sector = static_cast<std::int32_t>((u >> 15) & 0x1Fu) % 24 + 1;
    const std::int32_t arm = ((u >> 20) & 1u) ? 2 : 1;
    hits.push_back(SatElementHit{
        .arm = arm,
        .ring = ring,
        .sector = sector,
        .energy_gev = static_cast<float>(adc) * 0.05f,
        .r_cm = 10.0f + static_cast<float>(ring) * 3.0f,
        .phi_rad = (static_cast<float>(sector) - 0.5f) * (kTwoPi / 24.0f),
        .z_cm = (arm == 1) ? 232.0f : -232.0f,
    });
  }
  return hits;
}

std::vector<SatCluster> reconstruct_clusters(std::span<const SatElementHit> hits,
                                             float min_energy_gev) {
  // Aggregate per arm (1 and 2) for forward Bhabha/tagger energy reconstruction
  std::vector<SatCluster> out;
  for (std::int32_t arm = 1; arm <= 2; ++arm) {
    double sum_e = 0.0, sum_r = 0.0, sum_z = 0.0, sum_c = 0.0, sum_s = 0.0;
    std::int32_t count = 0;
    for (const auto& h : hits) {
      if (h.arm != arm || h.energy_gev <= 0.0f) continue;
      const double e = static_cast<double>(h.energy_gev);
      sum_e += e;
      sum_r += e * static_cast<double>(h.r_cm);
      sum_z += e * static_cast<double>(h.z_cm);
      sum_c += e * std::cos(static_cast<double>(h.phi_rad));
      sum_s += e * std::sin(static_cast<double>(h.phi_rad));
      ++count;
    }
    if (count > 0 && sum_e >= static_cast<double>(min_energy_gev)) {
      double phi = std::atan2(sum_s, sum_c);
      if (phi < 0.0) phi += 2.0 * std::numbers::pi_v<double>;
      out.push_back(SatCluster{
          .arm = arm,
          .energy_gev = static_cast<float>(sum_e),
          .r_cm = static_cast<float>(sum_r / sum_e),
          .phi_rad = static_cast<float>(phi),
          .z_cm = static_cast<float>(sum_z / sum_e),
          .num_hits = count,
      });
    }
  }
  std::sort(out.begin(), out.end(),
            [](const SatCluster& a, const SatCluster& b) {
              return a.energy_gev > b.energy_gev;
            });
  return out;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const SatCluster> clusters,
                                       std::int32_t detector_id) {
  std::vector<TanagraId> ids;
  ids.reserve(clusters.size());
  for (const auto& c : clusters) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = detector_id,
        .submod_id = c.arm,
        .descriptor = 125,
        .fwords = {c.r_cm, c.r_cm * c.phi_rad, c.z_cm, c.energy_gev},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::sat
