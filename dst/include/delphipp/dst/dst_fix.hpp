// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 post-DELANA Short-DST energy calibration, TPC ABST0 z-shift,
// and noise-cluster pruning (replacing simana_v94d/src/car/dstfix31.car
// & ecorrxx.car).

#pragma once

#include "delphipp/dst/short_dst.hpp"

namespace delphipp::dst {

struct DstFixConfig {
  float ecal_scale = 1.0f;
  float hcal_scale = 1.0f;
  float tpc_z_shift_cm = 0.0f;
  float min_neutral_energy_gev = 0.10f;
};

void apply_dst_fixes(DstEvent& ev, const DstFixConfig& cfg = {});

}  // namespace delphipp::dst
