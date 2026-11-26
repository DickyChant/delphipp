// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/vd.hpp"

#include <algorithm>
#include <cmath>
#include <tuple>

namespace delphipp::reco::vd {

std::vector<VdStripHit> unpack_raw_blocklet(std::span<const std::int32_t> words) {
  std::vector<VdStripHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::uint32_t adc = u & 0xFFu;
    if (adc == 0u) {
      continue;
    }
    hits.push_back(VdStripHit{
        .layer = static_cast<std::int32_t>((u >> 23) & 0x3u) % 3 + 1,
        .module = static_cast<std::int32_t>((u >> 18) & 0x1Fu) % 24 + 1,
        .strip = static_cast<std::int32_t>((u >> 8) & 0x3FFu),
        .pulse_height = static_cast<float>(adc),
        .noise = 2.0f,
        .is_z_side = ((u >> 25) & 1u) != 0u,
    });
  }
  return hits;
}

std::vector<VdCluster> find_clusters(std::span<const VdStripHit> strips,
                                     const VdConfig& cfg) {
  std::vector<VdStripHit> active;
  active.reserve(strips.size());
  for (const auto& s : strips) {
    if (s.noise > 0.0f && (s.pulse_height / s.noise) >= cfg.neighbor_sn_cut) {
      active.push_back(s);
    }
  }

  std::sort(active.begin(), active.end(), [](const VdStripHit& a, const VdStripHit& b) {
    return std::tie(a.layer, a.module, a.is_z_side, a.strip) <
           std::tie(b.layer, b.module, b.is_z_side, b.strip);
  });

  std::vector<VdCluster> clusters;
  std::size_t i = 0;
  while (i < active.size()) {
    std::size_t j = i + 1;
    while (j < active.size() &&
           active[j].layer == active[i].layer &&
           active[j].module == active[i].module &&
           active[j].is_z_side == active[i].is_z_side &&
           active[j].strip == active[j - 1].strip + 1) {
      ++j;
    }

    double sum_q = 0.0;
    double sum_qs = 0.0;
    double sum_n2 = 0.0;
    float max_sn = 0.0f;
    for (std::size_t k = i; k < j; ++k) {
      const float sn = active[k].pulse_height / active[k].noise;
      max_sn = std::max(max_sn, sn);
      sum_q += static_cast<double>(active[k].pulse_height);
      sum_qs += static_cast<double>(active[k].pulse_height) *
                static_cast<double>(active[k].strip);
      sum_n2 += static_cast<double>(active[k].noise) *
                static_cast<double>(active[k].noise);
    }

    if (max_sn >= cfg.seed_sn_cut && sum_q > 0.0) {
      const float centroid = static_cast<float>(sum_qs / sum_q);
      clusters.push_back(VdCluster{
          .layer = active[i].layer,
          .module = active[i].module,
          .is_z_side = active[i].is_z_side,
          .centroid_strip = centroid,
          .local_coord_cm = centroid * cfg.strip_pitch_cm,
          .total_charge = static_cast<float>(sum_q),
          .sn_ratio = static_cast<float>(sum_q / std::sqrt(sum_n2)),
          .width = static_cast<std::int32_t>(j - i),
      });
    }
    i = j;
  }

  return clusters;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const VdCluster> clusters) {
  std::vector<TanagraId> ids;
  ids.reserve(clusters.size());
  for (const auto& c : clusters) {
    const float r = (c.layer >= 1 && c.layer <= 3)
                        ? kLayerRadiiCm[static_cast<std::size_t>(c.layer - 1)]
                        : kLayerRadiiCm[0];
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = c.layer * 100 + c.module,
        .descriptor = c.is_z_side ? 126 : 125,
        .fwords = {r, c.local_coord_cm, c.total_charge, c.sn_ratio},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::vd
