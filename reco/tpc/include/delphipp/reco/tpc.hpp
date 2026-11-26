// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Time Projection Chamber (TPC) reconstruction
// (simana_v94d/src/car/tpcana.car): pad-row z-reconstruction using the exact
// 94D ABST0 calibration (1.10208 time bins vs 1.147 legacy, shifting z by
// +-232 um), plus 80% truncated-mean dE/dx computation.

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::tpc {

inline constexpr std::int32_t kDetectorId = 5;     // TPC detector code
inline constexpr float kEndplateZCm = 134.0f;      // Nominal TPC drift endplate |z| [cm]
inline constexpr float kTimeBinUs = 0.0735f;       // ~73.5 ns FADC time bin (3.3 ns / 0.04492 bins)
inline constexpr float kAbst0_94D = 1.10208f;      // 94D real-data TPC absolute T0 (ddb/tpct0_94c.asc)
inline constexpr float kAbst0_Sim = 1.147f;        // Simulation snapshot ABST0

struct TpcCalibration {
  float drift_velocity_cm_per_us = 6.72f;  // SLOW word 1
  float abst0_bins = kAbst0_94D;           // SLOW word 14
  float endplate_z_cm = kEndplateZCm;
  float time_bin_us = kTimeBinUs;
};

struct TpcPadHit {
  std::int32_t endplate = 0;  // 0 = ENP0 (+z), 1 = ENP1 (-z)
  std::int32_t sector = 1;    // 1..6 per endplate
  std::int32_t pad_row = 1;   // 1..16
  float r_cm = 50.0f;
  float r_phi_cm = 0.0f;
  float time_bin = 0.0f;      // Measured FADC drift time in bins
  float charge = 0.0f;
};

struct TpcSpacePoint {
  std::int32_t endplate = 0;
  std::int32_t sector = 1;
  std::int32_t pad_row = 1;
  float r_cm = 0.0f;
  float r_phi_cm = 0.0f;
  float z_cm = 0.0f;
  float charge = 0.0f;
};

// Unpack FADC raw-data blocklet words (tpcana.car) into TpcPadHit objects.
[[nodiscard]] std::vector<TpcPadHit> unpack_raw_blocklet(
    std::span<const std::int32_t> words,
    const TpcCalibration& calib = {});

// Reconstruct TPC z coordinate from FADC time bin and ABST0 calibration.
[[nodiscard]] float reconstruct_z(std::int32_t endplate,
                                  float time_bin,
                                  const TpcCalibration& calib) noexcept;

// Convert raw pad hits into 3D TPC space points.
[[nodiscard]] std::vector<TpcSpacePoint> reconstruct_space_points(
    std::span<const TpcPadHit> pad_hits,
    const TpcCalibration& calib);

// DELPHI TPC truncated-mean dE/dx (keeps lowest `keep_fraction`, default 80%).
[[nodiscard]] float truncated_mean_dedx(std::span<const float> wire_samples,
                                        float keep_fraction = 0.80f);

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const TpcSpacePoint> points);

}  // namespace delphipp::reco::tpc
