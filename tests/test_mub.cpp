// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>

#include "delphipp/reco/mub.hpp"

int main() {
  using namespace delphipp::reco::mub;

  // 1. Exhaustive round-trip over all 1,372 MUB chambers (EMCHAN <-> EMSLPB)
  for (std::int32_t idx = 1; idx <= kNumChambers; ++idx) {
    const auto loc = chamber_location(idx);
    assert(loc.has_value());
    const auto back = chamber_index(loc->sector, loc->layer, loc->pos);
    assert(back.has_value());
    assert(*back == idx);
  }
  // Invalid combinations
  assert(!chamber_index(49, 1, 1).has_value());  // leg region has no inner layers 1..3
  assert(!chamber_index(1, 5, 4).has_value());   // layer 5 max pos is 3

  // 2. EMCHID check for leg sector 10 -> sector 50 (offset 1344 + 14)
  assert(chamber_id(2, 10, 3, 1) == 1359);

  // 3. EMLTIM fine + coarse LTD timing
  LtdCalibration ltd{.coarse_step_ns = 8.0f, .fine_step_ns = {1.0f, 1.0f, 1.0f}};
  // raw = (10 << 2) | 0b11 = 43 -> ib0=1, ib1=1, jtime=40, tcors=320, tfine=6 -> 326 ns
  assert(std::fabs(ltd_time_ns(43, ltd) - 326.0f) < 1e-5f);
  // raw = (10 << 2) | 0b00 = 40 -> tcors=320, tfine=0 -> 320 ns
  assert(std::fabs(ltd_time_ns(40, ltd) - 320.0f) < 1e-5f);

  // 4. EMDFTT triplet reconstruction (real data with atmospheric pressure correction)
  MubChamberCalib ch{};
  ch.length_cm = 365.0f;
  ch.half_width_cm = 10.4f;
  ch.tzero_ns = {100.0f, 100.0f, 100.0f};
  ch.delay_poly = {0.0f, 0.5f, 0.0f, 0.0f, 0.0f};
  ch.anode_poly = {0.0f, 0.005f, 0.0f, 0.0f};

  MubConfig cfg{};
  cfg.vanode_cm_per_ns = 25.0f;
  cfg.veldrf_cm_per_ns = 0.005f;
  cfg.atmosp_mbar = 980.0f;  // 980/980 = 1.0

  // Symmetric near/far times -> zloc = 0 (center) -> z_cargo = 182.5 cm
  // tpr = 182.5 / 25.0 = 7.3 ns; td = 207.3 - 100 - 7.3 = 100 ns (>= 50)
  // ddrft = 0.117 + 0.005 * 100 = 0.617 cm
  const auto res = reconstruct_triplet(1, {207.3f, 500.0f, 500.0f}, ch, cfg);
  assert(res.has_value());
  assert(std::fabs(res->z_local_cm - 182.5f) < 1e-4f);
  assert(std::fabs(res->anode_time_ns - 100.0f) < 1e-4f);
  assert(std::fabs(res->drift_cm - 0.617f) < 1e-4f);

  // Reject unphysical time (near < anode)
  assert(!reconstruct_triplet(1, {500.0f, 400.0f, 500.0f}, ch, cfg).has_value());

  std::cout << "test_mub: PASS\n";
  return 0;
}
