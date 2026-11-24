// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delphipp/reco/muf.hpp"

int main() {
  using namespace delphipp::reco::muf;

  // 1. Time slew (EQTSLW): zero for simulation, non-zero cubic for real data
  MufConfig sim_cfg;
  sim_cfg.is_simulation = true;
  assert(time_slew(1, 450.0f, 0.2f, 1000.0f, sim_cfg) == 0.0f);

  MufConfig real_cfg;
  // At center of long chamber (x = 0), TSL == CTSLW0(1) = -7.680
  assert(std::fabs(time_slew(1, 400.0f, 0.2f, 1000.0f, real_cfg) - (-7.680f)) < 1e-4f);
  // At center of short chamber (x = 0), TSL == CTSLW0(2) = -5.12
  assert(std::fabs(time_slew(2, 300.0f, 0.2f, 750.0f, real_cfg) - (-5.12f)) < 1e-4f);

  // 2. Delay line distance (EQDLDF): 2-end symmetric hit at center
  MufRawTriplet h2{
      .quadrant = 1,
      .chamber = 4,
      .layer = 1,
      .ta_ns = 120.0f,
      .tl1_ns = 1100.0f,
      .tl2_ns = 1100.0f,
      .tca_ns = 10.0f,
      .tc1_ns = 0.0f,
      .tc2_ns = 0.0f,
      .vline_cm_per_ns = 0.2f,
      .chamber_len_cm = 400.0f,
      .time_sum_ns = 0.0f,
  };
  const auto dl2 = delay_line_distance(h2, 25.0f, real_cfg);
  assert(std::fabs(dl2.x_on_line_cm - 200.0f) < 1e-4f);
  assert(std::fabs(dl2.sigma_cm - 0.5f) < 1e-5f);

  // Single-end fallback scales error by 1.414
  MufRawTriplet h1 = h2;
  h1.tl2_ns = -1.0f;
  const auto dl1 = delay_line_distance(h1, 25.0f, real_cfg);
  assert(std::fabs(dl1.sigma_cm - 1.414f * 0.5f) < 1e-4f);

  // 3. Drift distance (EQDRIF): parabolic (< T1) and linear (>= T1)
  MufCalibration calib{};
  calib.vanode_cm_per_ns = 25.0f;
  calib.wire_pitch_cm = 18.0f;
  calib.qgeom[0][2] = 50.0f;   // x0
  calib.qgeom[0][4] = 0.0f;    // angle = 0
  calib.sx0par[0][0] = 1.0f;   // X0PARA = 1 cm
  calib.svdrif[0][0] = 0.005f; // VDRIFT = 0.005 cm/ns -> T1 = 100 ns

  // Parabolic regime: TDR = 25 ns = T1/4 -> DDR = X0 * sqrt(1/4) = 0.5 cm
  assert(std::fabs(drift_distance(1, 1, 25.0f, calib) - 0.5f) < 1e-5f);
  // Linear regime: TDR = 200 ns -> DDR = 1.0 + 0.005 * 100 = 1.5 cm
  assert(std::fabs(drift_distance(1, 1, 200.0f, calib) - 1.5f) < 1e-5f);

  // 4. Full hit reconstruction (EQSPAF)
  std::vector<MufRawTriplet> triplets = {h2};
  const auto hits = reconstruct_hits(triplets, calib, real_cfg);
  assert(hits.size() == 1);
  assert(hits[0].encoded_id == 1041);  // quad=1, chamber=4, layer=1
  // xwire = 50 + (4-1)*18 = 104 cm
  assert(std::fabs(hits[0].x_wire_cm - 104.0f) < 1e-4f);
  // tdr = 120 - 20 - (200/25) - 10 = 82 ns (< T1=100) -> parabolic
  const float expected_ddr = 1.0f * std::sqrt(82.0f / 100.0f);
  assert(std::fabs(hits[0].drift_dist_cm - expected_ddr) < 1e-4f);

  std::cout << "test_muf: PASS\n";
  return 0;
}
