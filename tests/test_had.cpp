// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

#include "delphipp/reco/had.hpp"
#include "delphipp/tanagra_store.hpp"

int main() {
  using namespace delphipp;
  using namespace delphipp::reco::had;

  // 1. ECALIB threshold & environmental corrections
  HadEnvironment mc_env{.run_number = 1000, .pressure_mbar = 990.0f,
                        .isobutane_pct = 30.5f, .temperature_c = 22.0f};
  assert(!calibrate_tower(1, 2, 500, mc_env).has_value());
  // MC (run <= 10581) ignores gas/pressure/temp: (4 + 0.5) * 500 * 0.001 = 2.25 GeV
  const auto e_mc = calibrate_tower(4, 2, 500, mc_env);
  assert(e_mc.has_value() && std::fabs(*e_mc - 2.25f) < 1e-5f);

  // Real data (run 47014 > 10581) applies pressure correction at 985 mbar
  HadEnvironment rd_env{.run_number = 47014, .pressure_mbar = 985.0f,
                        .isobutane_pct = 30.0f, .temperature_c = 20.0f};
  const auto e_rd = calibrate_tower(4, 2, 500, rd_env);
  const float expected_rd = 2.25f * (1.0f + 0.008056f * 10.0f);
  assert(e_rd.has_value() && std::fabs(*e_rd - expected_rd) < 1e-4f);

  // 2. Cyclic phi wraparound clustering (iphi=95 adjacent to iphi=0)
  const float two_pi = 2.0f * std::numbers::pi_v<float>;
  std::vector<HadTowerHit> towers = {
      {.itheta = 20, .iphi = 95, .ilayer = 1, .energy_gev = 3.0f,
       .r_cm = 350.0f, .theta_rad = 1.5f, .phi_rad = two_pi - 0.02f},
      {.itheta = 20, .iphi = 0, .ilayer = 2, .energy_gev = 3.0f,
       .r_cm = 370.0f, .theta_rad = 1.5f, .phi_rad = 0.02f},
  };
  const auto clusters = find_clusters(towers);
  assert(clusters.size() == 1);
  assert(clusters[0].num_towers == 2);
  assert(clusters[0].max_layer == 2);
  assert(std::fabs(clusters[0].energy_gev - 6.0f) < 1e-4f);
  // Circular mean of -0.02 and +0.02 is 0.0 (or 2*pi)
  assert(clusters[0].phi_rad < 1e-3f || std::fabs(clusters[0].phi_rad - two_pi) < 1e-3f);

  TanagraStore store;
  const auto ids = write_ter_banks(store, clusters);
  assert(ids.size() == 1);
  assert(store.list_ter(kDetectorId).size() == 1);

  std::cout << "test_had: PASS\n";
  return 0;
}
