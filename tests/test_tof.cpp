// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delphipp/reco/tof.hpp"
#include "delphipp/tanagra_store.hpp"

int main() {
  using namespace delphipp;
  using namespace delphipp::reco::tof;

  // 1. Sector poses: side A (z = +11.5) vs side C (z = -11.5)
  const auto s1 = sector_pose(1);
  assert(std::fabs(s1.z_cm - 11.5f) < 1e-5f);
  assert(s1.phi_rad >= 0.0f);
  const auto s30 = sector_pose(30);
  assert(std::fabs(s30.z_cm - (-11.5f)) < 1e-5f);

  // 2. EFSPA TD creation + pre-95 multihit sector deduplication
  TanagraStore store;
  std::vector<TofCounterHit> hits = {
      {.counter_num = 5, .r_cm = 311.25f, .r_phi_cm = 100.0f, .z_cm = 11.5f,
       .pm_ext = 12.0f, .pm_int = 11.5f},  // sector (5+3)/4 = 2
  };
  std::vector<std::int32_t> multihit = {2, 10};  // sector 2 already has hit -> skip; sector 10 added
  const auto ids = create_td_banks(store, hits, multihit, 94);
  assert(ids.size() == 2);
  assert(store.list_tdr(kDetectorId).size() == 2);
  assert(store.get_tdr(ids[0])->submod_id == 5);
  assert(store.get_tdr(ids[1])->submod_id == 192 + 100 * 10);

  // 3. EFNBUN majority-logic bunch identification (bunch 2 with bunint=247 ns -> [977, 1047])
  std::vector<float> tmaj2 = {1020.0f, 1005.0f, 1100.0f};
  assert(identify_bunch_majority(tmaj2, 4, 247.0f) == 2);
  assert(identify_bunch_majority({}, 4, 247.0f) == 0);

  std::cout << "test_tof: PASS\n";
  return 0;
}
