// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

#include "delphipp/reco/mus.hpp"

int main() {
  using namespace delphipp;
  using namespace delphipp::reco::mus;

  // 1. Amplitude shaping (EDAMSH): strip threshold 30, anode threshold 16
  assert(shape_amplitude(10, 30) < 0.0f);
  assert(std::fabs(shape_amplitude(10, 31) - 31.0f) < 1e-5f);
  assert(shape_amplitude(200, 16) < 0.0f);
  assert(std::fabs(shape_amplitude(200, 17) - 17.0f) < 1e-5f);

  // 2. Clustering (EDCLSR): 2 adjacent strips (10,11) + 1 isolated strip (20), 2 adjacent tubes (5,6)
  std::vector<MusRawHit> raw_hits = {
      {.det = 3, .layer = 1, .element = 10, .amplitude = 40.0f},
      {.det = 3, .layer = 1, .element = 11, .amplitude = 60.0f},
      {.det = 3, .layer = 1, .element = 20, .amplitude = 50.0f},
      {.det = 3, .layer = 2, .element = 5, .amplitude = 25.0f},
      {.det = 3, .layer = 2, .element = 6, .amplitude = 75.0f},
  };
  const auto clusters = find_clusters(raw_hits);
  assert(clusters.size() == 3);
  // Strip cluster 1: CM = (40*10 + 60*11)/100 = 10.6
  assert(clusters[0].layer == 1);
  assert(std::fabs(clusters[0].center_of_mass - 10.6f) < 1e-5f);
  assert(std::fabs(clusters[0].total_charge - 100.0f) < 1e-5f);
  // Strip cluster 2: CM = 20.0
  assert(std::fabs(clusters[1].center_of_mass - 20.0f) < 1e-5f);
  // Tube cluster 1: CM = (25*5 + 75*6)/100 = 5.75
  assert(clusters[2].layer == 2);
  assert(std::fabs(clusters[2].center_of_mass - 5.75f) < 1e-5f);

  // 3. Space point combination (EDCBCL) with 2 ambiguities
  MusGeometry geom{};
  geom.origin_ts1[2][0] = {0.0f, 0.0f, 0.0f};
  geom.pitch_ts1[2][0] = {4.0f, 0.0f};   // 4 cm strip pitch
  geom.origin_ts1[2][1] = {0.0f, 0.0f, 0.0f};
  geom.pitch_ts1[2][1] = {0.0f, 8.4f};   // 8.4 cm tube pitch
  // Identity rotation + translation to bottom module (x=0, y=-450, z=120) -> phi = 3pi/2
  geom.local_to_delphi[2][0] = {
      0.0f, -450.0f, 120.0f,
      1.0f, 0.0f, 0.0f,
      0.0f, 1.0f, 0.0f,
      0.0f, 0.0f, 1.0f,
  };

  const auto space_points = combine_clusters(clusters, geom);
  assert(space_points.size() == 1);
  assert(space_points[0].submod_id == 31);  // det=3, plane=1
  assert(space_points[0].solutions.size() == 2);
  // Solution 0: xloc = (10.6 - 1)*4 = 38.4, yloc = (5.75 - 1)*8.4 = 39.9, zloc = 0
  assert(std::fabs(space_points[0].solutions[0].local_xyz[0] - 38.4f) < 1e-4f);
  assert(std::fabs(space_points[0].solutions[0].local_xyz[1] - 39.9f) < 1e-4f);
  assert(std::fabs(space_points[0].solutions[0].local_xyz[2]) < 1e-5f);

  // 4. TE building (EDBLTE)
  const auto tes = build_te(space_points, geom);
  assert(tes.size() == 2);
  assert(tes[0].detector_id == kDetectorId);
  assert(tes[0].submod_id == 31);
  assert(std::fabs(tes[0].z - 120.0f) < 1e-4f);
  // Bottom module (y < 0, phi near 3pi/2) -> beta = 0.12 -> cov_z = (2.5 / cos(0.12))^2
  const float expected_cov_z = std::pow(2.5f / std::cos(0.12f), 2.0f);
  assert(std::fabs(tes[0].cov_z - expected_cov_z) < 1e-4f);

  std::cout << "test_mus: PASS\n";
  return 0;
}
