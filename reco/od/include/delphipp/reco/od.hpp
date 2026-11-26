// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Outer Detector (OD) reconstruction
// (simana_v94d/src/car/odana.car: 24 planks of 5 staggered drift-tube layers
// at R ~ 197..206 cm with fast z from dual-end timing).

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::od {

inline constexpr std::int32_t kDetectorId = 8;    // OD detector code
inline constexpr std::int32_t kNumPlanks = 24;
inline constexpr float kNominalRadiusCm = 202.0f;

struct OdTubeHit {
  std::int32_t plank = 1;  // 1..24
  std::int32_t layer = 1;  // 1..5
  float r_cm = kNominalRadiusCm;
  float r_phi_cm = 0.0f;
  float t_end_a_ns = 0.0f;
  float t_end_c_ns = 0.0f;
  bool has_z_timing = true;
};

struct OdTrackElement {
  std::int32_t plank = 1;
  float r_cm = kNominalRadiusCm;
  float r_phi_cm = 0.0f;
  float z_cm = 0.0f;
  std::int32_t num_layers = 0;
};

[[nodiscard]] std::vector<OdTubeHit> unpack_raw_blocklet(
    std::span<const std::int32_t> words);

[[nodiscard]] std::vector<OdTrackElement> reconstruct_planks(
    std::span<const OdTubeHit> hits,
    float v_prop_cm_per_ns = 20.0f);

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const OdTrackElement> elements);

}  // namespace delphipp::reco::od
