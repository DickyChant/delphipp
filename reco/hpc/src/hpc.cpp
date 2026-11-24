// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/hpc.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <queue>

namespace delphipp::reco::hpc {

namespace {

float delta_phi(float a, float b) noexcept {
  constexpr float kPi = std::numbers::pi_v<float>;
  float d = a - b;
  while (d > kPi) d -= 2.0f * kPi;
  while (d < -kPi) d += 2.0f * kPi;
  return d;
}

}  // namespace

std::vector<HpcShower> reconstruct_showers(std::span<const HpcChargeHit> hits,
                                           const HpcConfig& cfg) {
  std::vector<std::size_t> active;
  active.reserve(hits.size());
  for (std::size_t i = 0; i < hits.size(); ++i) {
    if (hits[i].energy_gev >= cfg.min_charge_gev) {
      active.push_back(i);
    }
  }

  std::vector<bool> visited(active.size(), false);
  std::vector<HpcShower> showers;

  for (std::size_t i = 0; i < active.size(); ++i) {
    if (visited[i]) continue;
    visited[i] = true;
    std::queue<std::size_t> q;
    q.push(i);
    std::vector<std::size_t> comp;

    while (!q.empty()) {
      const std::size_t cur = q.front();
      q.pop();
      comp.push_back(active[cur]);
      const auto& hc = hits[active[cur]];

      for (std::size_t j = 0; j < active.size(); ++j) {
        if (visited[j]) continue;
        const auto& hn = hits[active[j]];
        const float r_avg = 0.5f * (hc.r_cm + hn.r_cm);
        const float d_rphi = std::fabs(delta_phi(hc.phi_rad, hn.phi_rad)) * r_avg;
        const float d_z = std::fabs(hc.z_cm - hn.z_cm);
        if (d_rphi <= cfg.cluster_radius_rphi_cm && d_z <= cfg.cluster_radius_z_cm) {
          visited[j] = true;
          q.push(j);
        }
      }
    }

    double sum_e = 0.0;
    double sum_r = 0.0;
    double sum_z = 0.0;
    double sum_cos = 0.0;
    double sum_sin = 0.0;
    std::array<float, kNumLayers> lay_e{};

    for (const std::size_t idx : comp) {
      const auto& h = hits[idx];
      const double e = static_cast<double>(h.energy_gev);
      sum_e += e;
      sum_r += e * static_cast<double>(h.r_cm);
      sum_z += e * static_cast<double>(h.z_cm);
      sum_cos += e * std::cos(static_cast<double>(h.phi_rad));
      sum_sin += e * std::sin(static_cast<double>(h.phi_rad));
      if (h.layer >= 1 && h.layer <= static_cast<std::int32_t>(kNumLayers)) {
        lay_e[static_cast<std::size_t>(h.layer - 1)] += h.energy_gev;
      }
    }

    if (sum_e < static_cast<double>(cfg.min_shower_energy_gev)) continue;

    double phi = std::atan2(sum_sin, sum_cos);
    if (phi < 0.0) phi += 2.0 * std::numbers::pi_v<double>;

    showers.push_back(HpcShower{
        .module = hits[comp[0]].module,
        .energy_gev = static_cast<float>(sum_e),
        .r_cm = static_cast<float>(sum_r / sum_e),
        .phi_rad = static_cast<float>(phi),
        .z_cm = static_cast<float>(sum_z / sum_e),
        .layer_energy_gev = lay_e,
        .num_hits = static_cast<std::int32_t>(comp.size()),
    });
  }

  std::sort(showers.begin(), showers.end(),
            [](const HpcShower& a, const HpcShower& b) {
              return a.energy_gev > b.energy_gev;
            });
  return showers;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const HpcShower> showers) {
  std::vector<TanagraId> ids;
  ids.reserve(showers.size());
  for (const auto& s : showers) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = s.module,
        .descriptor = 125,
        .fwords = {s.r_cm, s.r_cm * s.phi_rad, s.z_cm, s.energy_gev},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::hpc
