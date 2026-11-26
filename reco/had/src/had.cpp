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

std::vector<HadTowerHit> unpack_raw_blocklet(std::span<const std::int32_t> words) {
  std::vector<HadTowerHit> hits;
  if (words.size() < 4) {
    return hits;
  }

  std::size_t istart = words.size();
  std::size_t nraw = 0;
  const std::uint32_t w2 = static_cast<std::uint32_t>(words[2]) & 0xFFFFu;
  if (w2 == 0x5A5Au || w2 == 0xA5A5u) {
    if (words[0] > 5) {
      nraw = static_cast<std::size_t>(words[0] - 5);
      istart = 4;
    }
  } else {
    for (std::size_t i = 1; i < words.size(); ++i) {
      const std::uint32_t lo = static_cast<std::uint32_t>(words[i]) & 0xFFFFu;
      const std::uint32_t hi = (static_cast<std::uint32_t>(words[i]) >> 16) & 0xFFFFu;
      if (lo == 0x5A5Au || lo == 0xA5A5u || hi == 0x5A5Au || hi == 0xA5A5u) {
        if (i >= 2 && words[i - 2] > 5) {
          istart = i + 2;
          nraw = static_cast<std::size_t>(words[i - 2] - 5);
          break;
        }
      }
    }
    if (istart >= words.size() && words[0] >= 5 &&
        static_cast<std::size_t>(words[0]) <= words.size()) {
      istart = 4;
      nraw = static_cast<std::size_t>(words[0] - 4);
    }
  }

  if (istart >= words.size()) {
    return hits;
  }
  const std::size_t end = std::min(words.size(), istart + nraw);
  constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
  for (std::size_t i = istart; i < end; ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::int32_t adc = static_cast<std::int32_t>(u & 0xFFu);
    if (adc <= 0) {
      continue;
    }
    const std::int32_t tower = static_cast<std::int32_t>((u >> 8) & 0xFFu);
    const std::int32_t port = static_cast<std::int32_t>((u >> 16) & 0x7Fu);
    const std::int32_t iphi = port % kNumPhiBins;
    const std::int32_t itheta = tower % 64;
    const std::int32_t ilayer = static_cast<std::int32_t>((u >> 24) & 0x3u) + 1;
    const float energy_gev = static_cast<float>(adc) * 0.17f;
    const float theta_rad = 0.25f + (static_cast<float>(itheta) + 0.5f) * (2.64f / 64.0f);
    const float phi_rad = (static_cast<float>(iphi) + 0.5f) * (kTwoPi / static_cast<float>(kNumPhiBins));
    const float r_cm = 320.0f + static_cast<float>(ilayer) * 15.0f;
    hits.push_back(HadTowerHit{
        .itheta = itheta,
        .iphi = iphi,
        .ilayer = ilayer,
        .energy_gev = energy_gev,
        .r_cm = r_cm,
        .theta_rad = theta_rad,
        .phi_rad = phi_rad,
    });
  }
  return hits;
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
