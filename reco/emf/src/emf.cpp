// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/emf.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <queue>

namespace delphipp::reco::emf {

std::vector<EmfBlockHit> unpack_raw_blocklet(std::span<const std::int32_t> words,
                                             const EmfConfig& cfg) {
  std::vector<EmfBlockHit> blocks;
  if (words.size() < 4) {
    return blocks;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::uint32_t adc = u & 0xFFFu;
    if (adc == 0u) {
      continue;
    }
    const float ph_Scale = (cfg.ph_per_gev_mc > 1.0f) ? cfg.ph_per_gev_mc : 1340.0f;
    const float energy_gev = static_cast<float>(adc) * 10.0f / ph_Scale;
    if (energy_gev < cfg.min_block_energy_gev || energy_gev > cfg.max_cluster_energy_gev) {
      continue;
    }
    const std::int32_t ix = static_cast<std::int32_t>((u >> 12) & 0x3Fu) - 32;
    const std::int32_t iy = static_cast<std::int32_t>((u >> 18) & 0x3Fu) - 32;
    const std::int32_t endcap = ((u >> 24) & 1u) ? 2 : 1;
    blocks.push_back(EmfBlockHit{
        .endcap = endcap,
        .ix = ix,
        .iy = iy,
        .energy_gev = energy_gev,
        .x_cm = static_cast<float>(ix) * kBlockPitchCm,
        .y_cm = static_cast<float>(iy) * kBlockPitchCm,
        .z_cm = (endcap == 1) ? kNominalZCm : -kNominalZCm,
    });
  }
  return blocks;
}

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
