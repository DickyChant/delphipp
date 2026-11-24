// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/had.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numbers>
#include <queue>

namespace delphipp::reco::had {

namespace {

std::int32_t cyclic_phi_dist(std::int32_t a, std::int32_t b) noexcept {
  const std::int32_t d = std::abs(a - b) % kNumPhiBins;
  return std::min(d, kNumPhiBins - d);
}

}  // namespace

std::optional<float> calibrate_tower(std::int32_t signal,
                                     std::int32_t threshold,
                                     std::int32_t tower_coeff_mev,
                                     const HadEnvironment& env) noexcept {
  if (signal < threshold) {
    return std::nullopt;
  }
  float e = static_cast<float>(signal) + 0.5f;
  if (e < 0.7f) {
    e = 0.85f;
  }
  float energy = e * static_cast<float>(tower_coeff_mev) * 0.001f * env.global_cal_factor;

  if (env.run_number > 10581) {
    if (env.pressure_mbar > 950.0f && env.pressure_mbar < 1000.0f) {
      energy *= (1.0f + 0.008056f * (env.pressure_mbar - 975.0f));
    }
    if (env.isobutane_pct > 29.0f && env.isobutane_pct < 31.0f) {
      energy *= (1.0f + 0.19f * (env.isobutane_pct - 30.0f));
    }
    if (env.temperature_c > 16.0f && env.temperature_c < 23.0f) {
      energy *= (1.0f - 0.030f * (env.temperature_c - 20.0f));
    }
  }
  return energy;
}

std::vector<HadCluster> find_clusters(std::span<const HadTowerHit> towers,
                                      float min_cluster_energy_gev) {
  std::vector<bool> visited(towers.size(), false);
  std::vector<HadCluster> clusters;

  for (std::size_t i = 0; i < towers.size(); ++i) {
    if (visited[i] || towers[i].energy_gev <= 0.0f) {
      continue;
    }
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
        if (visited[j] || towers[j].energy_gev <= 0.0f) {
          continue;
        }
        const auto& tn = towers[j];
        if (std::abs(tn.itheta - tc.itheta) <= 1 &&
            cyclic_phi_dist(tn.iphi, tc.iphi) <= 1 &&
            std::abs(tn.ilayer - tc.ilayer) <= 1) {
          visited[j] = true;
          q.push(j);
        }
      }
    }

    double sum_e = 0.0;
    double sum_r = 0.0;
    double sum_th = 0.0;
    double sum_cos = 0.0;
    double sum_sin = 0.0;
    std::int32_t max_lay = 0;
    for (const std::size_t idx : comp) {
      const auto& t = towers[idx];
      const double e = static_cast<double>(t.energy_gev);
      sum_e += e;
      sum_r += e * static_cast<double>(t.r_cm);
      sum_th += e * static_cast<double>(t.theta_rad);
      sum_cos += e * std::cos(static_cast<double>(t.phi_rad));
      sum_sin += e * std::sin(static_cast<double>(t.phi_rad));
      max_lay = std::max(max_lay, t.ilayer);
    }

    if (sum_e < static_cast<double>(min_cluster_energy_gev)) {
      continue;
    }

    double phi = std::atan2(sum_sin, sum_cos);
    if (phi < 0.0) {
      phi += 2.0 * std::numbers::pi_v<double>;
    }

    clusters.push_back(HadCluster{
        .energy_gev = static_cast<float>(sum_e),
        .r_cm = static_cast<float>(sum_r / sum_e),
        .theta_rad = static_cast<float>(sum_th / sum_e),
        .phi_rad = static_cast<float>(phi),
        .num_towers = static_cast<std::int32_t>(comp.size()),
        .max_layer = max_lay,
    });
  }

  std::sort(clusters.begin(), clusters.end(),
            [](const HadCluster& a, const HadCluster& b) {
              return a.energy_gev > b.energy_gev;
            });
  return clusters;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const HadCluster> clusters) {
  std::vector<TanagraId> ids;
  ids.reserve(clusters.size());
  for (const auto& c : clusters) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = c.max_layer,
        .descriptor = 126,  // Spherical (r, theta, phi) + energy
        .fwords = {c.r_cm, c.theta_rad, c.phi_rad, c.energy_gev},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::had
