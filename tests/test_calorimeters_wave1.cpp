// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

#include "delphipp/reco/hpc.hpp"
#include "delphipp/reco/sat.hpp"
#include "delphipp/reco/stc.hpp"
#include "delphipp/tanagra_store.hpp"

int main() {
  using namespace delphipp;

  TanagraStore store;

  // 1. HPC barrel ECAL shower reconstruction
  std::vector<reco::hpc::HpcChargeHit> hpc_hits = {
      {.module = 5, .layer = 2, .r_cm = 215.0f, .phi_rad = 1.00f, .z_cm = 30.0f, .energy_gev = 12.0f},
      {.module = 5, .layer = 4, .r_cm = 225.0f, .phi_rad = 1.01f, .z_cm = 31.0f, .energy_gev = 30.0f},
      {.module = 5, .layer = 6, .r_cm = 235.0f, .phi_rad = 1.00f, .z_cm = 30.5f, .energy_gev = 3.6f},
  };
  const auto hpc_showers = reco::hpc::reconstruct_showers(hpc_hits);
  assert(hpc_showers.size() == 1);
  assert(std::fabs(hpc_showers[0].energy_gev - 45.6f) < 1e-3f);
  assert(std::fabs(hpc_showers[0].layer_energy_gev[3] - 30.0f) < 1e-4f);
  assert(reco::hpc::write_ter_banks(store, hpc_showers).size() == 1);

  // 2. STIC luminosity Bhabha candidate (back-to-back 45 GeV arms)
  const float pi = std::numbers::pi_v<float>;
  std::vector<reco::stc::StcTowerHit> stc_hits = {
      {.arm = 1, .ring = 5, .sector = 1, .energy_gev = 44.5f, .r_cm = 18.0f, .phi_rad = 0.2f, .z_cm = 220.0f},
      {.arm = 2, .ring = 5, .sector = 9, .energy_gev = 45.1f, .r_cm = 18.0f, .phi_rad = pi + 0.25f, .z_cm = -220.0f},
  };
  const auto stc_clus = reco::stc::find_clusters(stc_hits);
  assert(stc_clus.size() == 2);
  assert(reco::stc::is_bhabha_candidate(stc_clus));
  assert(reco::stc::write_ter_banks(store, stc_clus).size() == 2);

  // 3. SAT forward tagger
  std::vector<reco::sat::SatElementHit> sat_hits = {
      {.arm = 1, .ring = 3, .sector = 4, .energy_gev = 15.0f, .r_cm = 22.0f, .phi_rad = 0.5f, .z_cm = 232.0f},
      {.arm = 1, .ring = 4, .sector = 4, .energy_gev = 10.0f, .r_cm = 24.0f, .phi_rad = 0.5f, .z_cm = 232.0f},
  };
  const auto sat_clus = reco::sat::reconstruct_clusters(sat_hits);
  assert(sat_clus.size() == 1);
  assert(std::fabs(sat_clus[0].energy_gev - 25.0f) < 1e-4f);
  assert(reco::sat::write_ter_banks(store, sat_clus).size() == 1);

  assert(store.ter_count() == 4);
  std::cout << "test_calorimeters_wave1: PASS\n";
  return 0;
}
