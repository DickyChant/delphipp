// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numbers>

#include "delphipp/reco/hof.hpp"
#include "delphipp/zebra_store.hpp"

namespace {

std::int32_t pack_bits(std::uint32_t val, int pos_0based) {
  return static_cast<std::int32_t>(val << pos_0based);
}

}  // namespace

int main() {
  using namespace delphipp;
  using namespace delphipp::reco::hof;

  OwnedZebraStore store(2048);
  auto view = store.view();

  // Populate calibration for Endcap 1, Quadrant 2, Counter 5 (1-based)
  HofCalibration calib;
  const std::size_t geo_idx = ((1 - 1) * 112 + (2 - 1) * 28 + (5 - 1)) * 5;
  calib.geoddb[geo_idx + 0] = 100.0f;   // x
  calib.geoddb[geo_idx + 1] = -100.0f;  // y (negative -> phi in [3pi/2, 2pi))
  calib.geoddb[geo_idx + 2] = 450.0f;   // z
  calib.geoddb[geo_idx + 3] = 6.0f;     // dx
  calib.geoddb[geo_idx + 4] = 9.0f;     // dy

  // Map raw channel bit 3 in Endcap 1, Quad 2 -> Counter 5
  const std::size_t pad_idx = (1 - 1) * 128 + (2 - 1) * 32 + 3;
  calib.ichan[pad_idx] = 5;

  // Map timing slot=7, chan=12 for Endcap 1, Quad 2
  const std::size_t time_base = 256 + 11 * (1 - 1);
  calib.ichan[time_base + 2] = 7;       // ISLDB
  calib.ichan[time_base + 2 + 5] = 12;  // ICHDB

  // Synthetic EGRAW bank at L = 100 (module=29, numcap=1, bit 3 set in quad 2)
  const std::int32_t lraw = 100;
  view.set_lq(lraw, 0);
  view.set_iq(lraw - 1, 10);  // nd > 5
  view.set_iq(lraw + 2, pack_bits(29u, 8) | pack_bits(1u, 16));
  view.set_iq(lraw + 4 + 2, pack_bits(1u, 3));  // Quad 2, bit 3

  // Synthetic EQRAW bank at L = 300 (module=30, endcap=1, 1 timing hit = 310)
  const std::int32_t mraw = 300;
  view.set_lq(mraw, 0);
  view.set_iq(mraw - 1, 40);
  view.set_iq(mraw + 2, pack_bits(30u, 8) | pack_bits(1u, 16));
  // Header word at mraw + 1 + 4: iwd1=10 (bits 0..16), nh=1 (bits 17..26), idrop=0
  view.set_iq(mraw + 5, pack_bits(10u, 0) | pack_bits(1u, 17));
  // Slot word at mraw + 10 + 3: islot=7 (1-based bit 17 -> 0-based bit 16)
  view.set_iq(mraw + 13, pack_bits(7u, 16));
  // Hit word at mraw + 10 + 5 + 1: ihoft=310 (bits 0..13), nchan=12 (1-based bit 15 -> 0-based bit 14)
  view.set_iq(mraw + 16, pack_bits(310u, 0) | pack_bits(12u, 14));

  HofStats stats;
  stats.reset_run();
  stats.begin_event();

  const HofTdResult td = unpack_td(view, lraw, mraw, calib, stats);
  assert(td.counter_hits == 1);
  assert(td.hits.size() == 2);  // 1 counter hit + 1 timing hit
  assert(stats.nsgqua[1] == 1);
  assert(stats.nsgful[1] == 1);
  assert(std::fabs(td.hits[0].id_word() - 1205.0f) < 1e-4f);
  assert(std::fabs(td.hits[1].id_word() - (310.0f * 100.0f + 12.0f)) < 1e-3f);

  const auto tes = build_te(view, mraw, td, calib);
  assert(tes.size() == 1);
  assert(tes[0].get_int(1) == kDetectorId);
  assert(tes[0].get_int(2) == 12);  // 1205 / 100
  assert(std::fabs(tes[0].get_float(10) - 100.0f) < 1e-5f);
  assert(std::fabs(tes[0].get_float(11) - (-100.0f)) < 1e-5f);
  assert(std::fabs(tes[0].get_float(12) - 450.0f) < 1e-5f);
  // phi = atan2(-100, 100) + 2*pi = 7*pi/4
  const float expected_phi = 1.75f * std::numbers::pi_v<float>;
  assert(std::fabs(tes[0].get_float(14) - expected_phi) < 1e-5f);
  assert(std::fabs(tes[0].get_float(16) - 310.0f) < 1e-4f);
  // Covariance: dx^2/3 = 12, dy^2/3 = 27
  assert(std::fabs(tes[0].get_float(17) - 12.0f) < 1e-5f);
  assert(std::fabs(tes[0].get_float(19) - 27.0f) < 1e-5f);

  // Bunch identification: 310 is within 24 counts of TIMM=307 -> bunch 1
  assert(identify_bunch(tes, 4, 247.0f) == 1);

  std::cout << "test_hof: PASS\n";
  return 0;
}
