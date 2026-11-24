// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Time-of-Flight (TOF) reconstruction
// (simana_v94d/src/car/tofana.car: EFSPA, EFNBUN).

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::tof {

inline constexpr std::int32_t kDetectorId = 13;    // ILITOF
inline constexpr float kTofRadiusCm = 311.25f;     // Nominal barrel TOF radius
inline constexpr std::int32_t kNumSectors = 48;    // 24 side A + 24 side C

struct TofCounterHit {
  std::int32_t counter_num = 1;  // 1..192 (NCTNUM)
  float r_cm = kTofRadiusCm;     // DINTER(i,2)
  float r_phi_cm = 0.0f;         // DINTER(i,3)
  float z_cm = 0.0f;             // DINTER(i,1)
  float pm_ext = 0.0f;           // DATPME(i)
  float pm_int = 0.0f;           // DATPMI(i)
};

// Sector azimuth [0, 2*pi) and nominal z (+11.5 for 1..24, -11.5 for 25..48)
struct TofSectorPose {
  float r_cm = kTofRadiusCm;
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
};

[[nodiscard]] TofSectorPose sector_pose(std::int32_t sector_1based) noexcept;

// Port of SUBROUTINE EFSPA: creates TOF TDR banks directly in TanagraStore
// for calibrated counter hits and (optionally) pre-95 multihit sectors.
// Returns the created TDR IDs.
std::vector<TanagraId> create_td_banks(
    TanagraStore& store,
    std::span<const TofCounterHit> hits,
    std::span<const std::int32_t> multihit_sectors = {},
    std::int32_t year_version = 94);

// Port of SUBROUTINE EFNBUN (majority-logic branch): returns 1-based bunch
// number (1..nbulep) or 0 if no window matches.
[[nodiscard]] std::int32_t identify_bunch_majority(
    std::span<const float> tmaj2_times_ns,
    std::int32_t nbulep,
    float bunint_ns) noexcept;

}  // namespace delphipp::reco::tof
