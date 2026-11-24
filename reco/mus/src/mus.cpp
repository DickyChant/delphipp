// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/mus.hpp"

#include <cmath>
#include <cstdlib>
#include <numbers>

namespace delphipp::reco::mus {

namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;

}  // namespace

float shape_amplitude(std::int32_t iadr,
                      std::int32_t iadc,
                      const MusConfig& cfg) noexcept {
  if (iadr < 192) {
    if (iadc <= cfg.strip_pedestal) {
      return -1.0f;
    }
  } else {
    if (iadc <= cfg.anode_pedestal) {
      return -1.0f;
    }
  }
  return static_cast<float>(iadc);
}

std::vector<MusCluster> find_clusters(std::span<const MusRawHit> raw_hits) {
  // Bucket raw hit indices by (det 0..15, layer 0..3)
  std::array<std::array<std::vector<std::size_t>, kNumLayers>, kNumDetectors> by_layer;
  for (std::size_t idx = 0; idx < raw_hits.size(); ++idx) {
    const auto& h = raw_hits[idx];
    if (h.det < 1 || h.det > static_cast<std::int32_t>(kNumDetectors) ||
        h.layer < 1 || h.layer > static_cast<std::int32_t>(kNumLayers)) {
      continue;
    }
    auto& bucket = by_layer[static_cast<std::size_t>(h.det - 1)]
                           [static_cast<std::size_t>(h.layer - 1)];
    if (bucket.size() < kMaxHitsPerLayer) {
      bucket.push_back(idx);
    }
  }

  std::vector<MusCluster> clusters;

  for (std::size_t d = 0; d < kNumDetectors; ++d) {
    for (std::size_t l = 0; l < kNumLayers; ++l) {
      const auto& bucket = by_layer[d][l];
      if (bucket.empty()) {
        continue;
      }

      std::size_t layer_cluster_count = 0;
      std::int32_t prev_elem = 0;
      float weighted_sum = 0.0f;
      float amp_sum = 0.0f;
      std::vector<std::size_t> current_hits;

      auto flush_cluster = [&]() {
        if (current_hits.empty() || amp_sum <= 0.0f) {
          return;
        }
        MusCluster clu;
        clu.det = static_cast<std::int32_t>(d + 1);
        clu.layer = static_cast<std::int32_t>(l + 1);
        clu.center_of_mass = weighted_sum / amp_sum;
        clu.total_charge = amp_sum;
        clu.raw_hit_indices = std::move(current_hits);
        clusters.push_back(std::move(clu));
        ++layer_cluster_count;
        weighted_sum = 0.0f;
        amp_sum = 0.0f;
        current_hits.clear();
      };

      for (std::size_t i = 0; i < bucket.size(); ++i) {
        const std::size_t raw_idx = bucket[i];
        const auto& h = raw_hits[raw_idx];
        if (i > 0 && std::abs(h.element - prev_elem) > 1) {
          flush_cluster();
          if (layer_cluster_count >= kMaxClustersPerLayer) {
            break;
          }
        }
        weighted_sum += h.amplitude * static_cast<float>(h.element);
        amp_sum += h.amplitude;
        if (current_hits.size() < kMaxHitsPerCluster) {
          current_hits.push_back(raw_idx);
        }
        prev_elem = h.element;
      }
      if (layer_cluster_count < kMaxClustersPerLayer) {
        flush_cluster();
      }
    }
  }

  return clusters;
}

std::vector<MusSpacePoint> combine_clusters(std::span<const MusCluster> clusters,
                                            const MusGeometry& geom) {
  std::array<std::array<std::vector<std::size_t>, kNumLayers>, kNumDetectors> clu_by_layer;
  for (std::size_t idx = 0; idx < clusters.size(); ++idx) {
    const auto& c = clusters[idx];
    if (c.det < 1 || c.det > static_cast<std::int32_t>(kNumDetectors) ||
        c.layer < 1 || c.layer > static_cast<std::int32_t>(kNumLayers)) {
      continue;
    }
    auto& vec = clu_by_layer[static_cast<std::size_t>(c.det - 1)]
                            [static_cast<std::size_t>(c.layer - 1)];
    if (vec.size() < kMaxClustersPerLayer) {
      vec.push_back(idx);
    }
  }

  std::vector<MusSpacePoint> points;

  for (std::size_t d = 0; d < kNumDetectors; ++d) {
    // Tube layers are 1-based 2 and 3 (0-based 1 and 2); complementary strip layers are 1 and 4 (0 and 3)
    for (std::size_t tube_layer = 1; tube_layer <= 2; ++tube_layer) {
      const std::size_t strip_layer = (tube_layer == 1) ? 0 : 3;
      const std::int32_t plane = (tube_layer == 1) ? 1 : 2;

      const auto& tube_clus = clu_by_layer[d][tube_layer];
      const auto& strip_clus = clu_by_layer[d][strip_layer];
      if (tube_clus.empty() || strip_clus.empty()) {
        continue;
      }

      for (const std::size_t tc_idx : tube_clus) {
        if (points.size() >= kMaxSpacePoints) {
          return points;
        }
        MusSpacePoint sp;
        sp.det = static_cast<std::int32_t>(d + 1);
        sp.plane = plane;
        sp.submod_id = sp.det * 10 + plane;

        const auto& tc = clusters[tc_idx];
        for (const std::size_t sc_idx : strip_clus) {
          if (sp.solutions.size() >= kMaxSolutionsPerPoint) {
            break;
          }
          const auto& sc = clusters[sc_idx];
          const float xloc = geom.origin_ts1[d][strip_layer][0] +
                             (sc.center_of_mass - 1.0f) * geom.pitch_ts1[d][strip_layer][0];
          const float yloc = geom.origin_ts1[d][tube_layer][1] +
                             (tc.center_of_mass - 1.0f) * geom.pitch_ts1[d][tube_layer][1];
          // musana.car line 1481: zloc = 0 (center of module plane)
          sp.solutions.push_back(MusSpacePointSolution{
              .local_xyz = {xloc, yloc, 0.0f},
              .tube_cluster_idx = tc_idx,
              .strip_cluster_idx = sc_idx,
          });
        }
        points.push_back(std::move(sp));
      }
    }
  }

  return points;
}

std::vector<MusTeBank> build_te(std::span<const MusSpacePoint> space_points,
                                const MusGeometry& geom,
                                const MusConfig& cfg) {
  std::vector<MusTeBank> tes;
  for (const auto& sp : space_points) {
    if (sp.det < 1 || sp.det > static_cast<std::int32_t>(kNumDetectors) ||
        sp.plane < 1 || sp.plane > 2) {
      continue;
    }
    const auto& rmatr = geom.local_to_delphi[static_cast<std::size_t>(sp.det - 1)]
                                            [static_cast<std::size_t>(sp.plane - 1)];

    for (const auto& sol : sp.solutions) {
      const Vec3 xyz = dtrans(sol.local_xyz, rmatr, 0);
      const float r = std::hypot(xyz[0], xyz[1]);
      float phi = 0.0f;
      if (xyz[0] != 0.0f) {
        phi = std::atan2(xyz[1], xyz[0]);
        if (phi < 0.0f) {
          phi += kTwoPi;
        }
      } else if (xyz[1] > 0.0f) {
        phi = kPi * 0.5f;
      } else if (xyz[1] < 0.0f) {
        phi = 1.5f * kPi;
      }

      const float phi1 = (std::fabs(xyz[0]) != 0.0f)
                             ? std::atan2(std::fabs(xyz[1]), std::fabs(xyz[0]))
                             : (kPi * 0.5f);
      const float alpha = (phi1 > kPi * 0.25f) ? (kPi * 0.5f - phi1) : phi1;
      const float beta = (phi >= 1.25f * kPi && phi <= 1.75f * kPi)
                             ? cfg.beta_bottom_rad
                             : cfg.beta_top_side_rad;

      const float err_rphi = cfg.sigma_strip_cm * std::cos(alpha);
      const float err_z = cfg.sigma_tube_cm / std::cos(beta);

      MusTeBank te;
      te.detector_id = kDetectorId;
      te.submod_id = sp.submod_id;
      te.descriptor = 125;
      te.r = r;
      te.r_phi = r * phi;
      te.z = xyz[2];
      te.cov_rphi = err_rphi * err_rphi;
      te.cov_rphi_z = 0.0f;
      te.cov_z = err_z * err_z;
      tes.push_back(te);
    }
  }
  return tes;
}

}  // namespace delphipp::reco::mus
