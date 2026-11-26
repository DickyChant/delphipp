// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/stc.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numbers>
#include <queue>

namespace delphipp::reco::stc {

namespace {

std::int32_t cyclic_sector_dist(std::int32_t a, std::int32_t b) noexcept {
  const std::int32_t d = std::abs(a - b) % kNumSectors;
  return std::min(d, kNumSectors - d);
}

}  // namespace

std::vector<StcTowerHit> unpack_raw_blocklet(std::span<const std::int32_t> words) {
  std::vector<StcTowerHit> towers;
  if (words.size() < 4) {
    return towers;
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
    const std::int32_t ring = static_cast<std::int32_t>((u >> 12) & 0xFu) % kNumRings + 1;
    const std::int32_t sector = static_cast<std::int32_t>((u >> 16) & 0xFu) % kNumSectors + 1;
    const std::int32_t arm = ((u >> 20) & 1u) ? 2 : 1;
    towers.push_back(StcTowerHit{
        .arm = arm,
        .ring = ring,
        .sector = sector,
        .energy_gev = static_cast<float>(adc) * 0.05f,
        .r_cm = 6.5f + static_cast<float>(ring) * 3.5f,
        .phi_rad = (static_cast<float>(sector) - 0.5f) * (kTwoPi / static_cast<float>(kNumSectors)),
        .z_cm = (arm == 1) ? kNominalZCm : -kNominalZCm,
    });
  }
  return towers;
}

std::vector<StcCluster> find_clusters(std::span<const StcTowerHit> towers,
                                      float min_cluster_energy_gev) {
  std::vector<bool> visited(towers.size(), false);
  std::vector<StcCluster> clusters;

  for (std::size_t i = 0; i < towers.size(); ++i) {
    if (visited[i] || towers[i].energy_gev <= 0.0f) continue;
    visited[i] = true;
    std::queue<std::size_t> q;
    q.push(i);
    std::vector<std::size_t> comp;

    while (!q.empty()) {
      const std::size_t cur = q.front();
      q.pop();
      comp.push_back(cur);
      const auto& tc = towers[cur];
      for (std::size_t j = 0; j < towers.size(); ++j) {
        if (visited[j] || towers[j].energy_gev <= 0.0f) continue;
        const auto& tn = towers[j];
        if (tn.arm == tc.arm &&
            std::abs(tn.ring - tc.ring) <= 1 &&
            cyclic_sector_dist(tn.sector, tc.sector) <= 1) {
          visited[j] = true;
          q.push(j);
        }
      }
    }

    double sum_e = 0.0, sum_r = 0.0, sum_z = 0.0, sum_c = 0.0, sum_s = 0.0;
    for (const std::size_t idx : comp) {
      const auto& t = towers[idx];
      const double e = static_cast<double>(t.energy_gev);
      sum_e += e;
      sum_r += e * static_cast<double>(t.r_cm);
      sum_z += e * static_cast<double>(t.z_cm);
      sum_c += e * std::cos(static_cast<double>(t.phi_rad));
      sum_s += e * std::sin(static_cast<double>(t.phi_rad));
    }
    if (sum_e < static_cast<double>(min_cluster_energy_gev)) continue;

    double phi = std::atan2(sum_s, sum_c);
    if (phi < 0.0) phi += 2.0 * std::numbers::pi_v<double>;

    clusters.push_back(StcCluster{
        .arm = towers[comp[0]].arm,
        .energy_gev = static_cast<float>(sum_e),
        .r_cm = static_cast<float>(sum_r / sum_e),
        .phi_rad = static_cast<float>(phi),
        .z_cm = static_cast<float>(sum_z / sum_e),
        .num_towers = static_cast<std::int32_t>(comp.size()),
    });
  }

  std::sort(clusters.begin(), clusters.end(),
            [](const StcCluster& a, const StcCluster& b) {
              return a.energy_gev > b.energy_gev;
            });
  return clusters;
}

bool is_bhabha_candidate(std::span<const StcCluster> clusters,
                         float min_arm_energy_gev,
                         float max_acoplanarity_rad) noexcept {
  const StcCluster* arm1 = nullptr;
  const StcCluster* arm2 = nullptr;
  for (const auto& c : clusters) {
    if (c.arm == 1 && (!arm1 || c.energy_gev > arm1->energy_gev)) arm1 = &c;
    if (c.arm == 2 && (!arm2 || c.energy_gev > arm2->energy_gev)) arm2 = &c;
  }
  if (!arm1 || !arm2) return false;
  if (arm1->energy_gev < min_arm_energy_gev || arm2->energy_gev < min_arm_energy_gev) {
    return false;
  }
  constexpr float kPi = std::numbers::pi_v<float>;
  float dphi = std::fabs(arm1->phi_rad - arm2->phi_rad);
  if (dphi > kPi) dphi = 2.0f * kPi - dphi;
  return std::fabs(kPi - dphi) <= max_acoplanarity_rad;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const StcCluster> clusters) {
  std::vector<TanagraId> ids;
  ids.reserve(clusters.size());
  for (const auto& c : clusters) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = c.arm,
        .descriptor = 125,
        .fwords = {c.r_cm, c.r_cm * c.phi_rad, c.z_cm, c.energy_gev},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::stc
