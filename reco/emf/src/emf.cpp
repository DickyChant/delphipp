// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/emf.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <queue>

namespace delphipp::reco::emf {

std::vector<EmfCluster> find_clusters(std::span<const EmfBlockHit> blocks,
                                      const EmfConfig& cfg) {
  std::vector<std::size_t> active;
  active.reserve(blocks.size());
  for (std::size_t i = 0; i < blocks.size(); ++i) {
    if (blocks[i].energy_gev >= cfg.min_block_energy_gev) {
      active.push_back(i);
    }
  }

  std::vector<bool> visited(active.size(), false);
  std::vector<EmfCluster> clusters;

  for (std::size_t i = 0; i < active.size(); ++i) {
    if (visited[i]) {
      continue;
    }
    visited[i] = true;
    std::queue<std::size_t> q;
    q.push(i);

    std::vector<std::size_t> comp;
    while (!q.empty()) {
      const std::size_t cur = q.front();
      q.pop();
      comp.push_back(active[cur]);
      const auto& bc = blocks[active[cur]];

      for (std::size_t j = 0; j < active.size(); ++j) {
        if (visited[j]) {
          continue;
        }
        const auto& bn = blocks[active[j]];
        if (bn.endcap == bc.endcap &&
            std::abs(bn.ix - bc.ix) <= 1 &&
            std::abs(bn.iy - bc.iy) <= 1) {
          visited[j] = true;
          q.push(j);
        }
      }
    }

    double sum_e = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    for (const std::size_t idx : comp) {
      const auto& b = blocks[idx];
      const double e = static_cast<double>(b.energy_gev);
      sum_e += e;
      sum_x += e * static_cast<double>(b.x_cm);
      sum_y += e * static_cast<double>(b.y_cm);
      sum_z += e * static_cast<double>(b.z_cm);
    }

    if (sum_e < static_cast<double>(cfg.min_cluster_energy_gev) ||
        sum_e > static_cast<double>(cfg.max_cluster_energy_gev)) {
      continue;
    }

    const double mean_x = sum_x / sum_e;
    const double mean_y = sum_y / sum_e;
    const double mean_z = sum_z / sum_e;

    double var_x = 0.0;
    double var_y = 0.0;
    for (const std::size_t idx : comp) {
      const auto& b = blocks[idx];
      const double w = static_cast<double>(b.energy_gev) / sum_e;
      const double dx = static_cast<double>(b.x_cm) - mean_x;
      const double dy = static_cast<double>(b.y_cm) - mean_y;
      var_x += w * dx * dx;
      var_y += w * dy * dy;
    }

    // Floor dispersion at single-block uniform variance (pitch / sqrt(12))
    const float min_sigma = kBlockPitchCm / std::sqrt(12.0f);
    clusters.push_back(EmfCluster{
        .endcap = blocks[comp[0]].endcap,
        .energy_gev = static_cast<float>(sum_e),
        .x_cm = static_cast<float>(mean_x),
        .y_cm = static_cast<float>(mean_y),
        .z_cm = static_cast<float>(mean_z),
        .sigma_x_cm = std::max(min_sigma, static_cast<float>(std::sqrt(var_x))),
        .sigma_y_cm = std::max(min_sigma, static_cast<float>(std::sqrt(var_y))),
        .num_blocks = static_cast<std::int32_t>(comp.size()),
    });
  }

  // Deterministic descending energy ordering
  std::sort(clusters.begin(), clusters.end(),
            [](const EmfCluster& a, const EmfCluster& b) {
              return a.energy_gev > b.energy_gev;
            });
  return clusters;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const EmfCluster> clusters) {
  std::vector<TanagraId> ids;
  ids.reserve(clusters.size());
  for (const auto& c : clusters) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = c.endcap,
        .descriptor = 124,  // Cartesian (x, y, z) + energy
        .fwords = {c.x_cm, c.y_cm, c.z_cm, c.energy_gev,
                   c.sigma_x_cm, c.sigma_y_cm},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::emf
