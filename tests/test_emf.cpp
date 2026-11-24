// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delphipp/reco/emf.hpp"
#include "delphipp/tanagra_store.hpp"

int main() {
  using namespace delphipp;
  using namespace delphipp::reco::emf;

  std::vector<EmfBlockHit> blocks = {
      // Endcap 1 cluster: 2 adjacent blocks (total E = 10 GeV, centroid x = 104.2)
      {.endcap = 1, .ix = 10, .iy = 10, .energy_gev = 2.5f, .x_cm = 100.0f, .y_cm = 50.0f, .z_cm = 284.0f},
      {.endcap = 1, .ix = 11, .iy = 10, .energy_gev = 7.5f, .x_cm = 105.6f, .y_cm = 50.0f, .z_cm = 284.0f},
      // Endcap 2 cluster at same (ix, iy) must NOT merge across endcaps
      {.endcap = 2, .ix = 10, .iy = 10, .energy_gev = 20.0f, .x_cm = 100.0f, .y_cm = 50.0f, .z_cm = -284.0f},
      // Below-threshold noise block
      {.endcap = 1, .ix = 20, .iy = 20, .energy_gev = 0.01f, .x_cm = 150.0f, .y_cm = 150.0f, .z_cm = 284.0f},
  };

  const auto clusters = find_clusters(blocks);
  assert(clusters.size() == 2);
  // Sorted descending by energy: Endcap 2 (20 GeV) first, Endcap 1 (10 GeV) second
  assert(clusters[0].endcap == 2);
  assert(std::fabs(clusters[0].energy_gev - 20.0f) < 1e-4f);
  assert(clusters[1].endcap == 1);
  assert(clusters[1].num_blocks == 2);
  assert(std::fabs(clusters[1].energy_gev - 10.0f) < 1e-4f);
  assert(std::fabs(clusters[1].x_cm - 104.2f) < 1e-4f);

  TanagraStore store;
  const auto ids = write_ter_banks(store, clusters);
  assert(ids.size() == 2);
  assert(store.list_ter(kDetectorId).size() == 2);

  std::cout << "test_emf: PASS\n";
  return 0;
}
