// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI ID Straw Tube Detector (STW / 1995+ upgrade)
// reconstruction (simana_v94d/src/car/stwana.car: 5 concentric 192-straw
// layers grouped into 12 azimuthal sectors, fitted into local track elements).

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::stw {

inline constexpr std::int32_t kDetectorId = 31;  // ID Straw detector code
inline constexpr std::int32_t kNumLayers = 5;
inline constexpr std::int32_t kNumStrawsPerLayer = 192;
inline constexpr std::int32_t kNumSectors = 12;
inline constexpr std::array<float, 5> kLayerRadiiCm = {23.2f, 24.4f, 25.6f, 26.8f, 28.0f};

struct StwHit {
  std::int32_t layer = 1;    // 1..5
  std::int32_t straw = 1;    // 1..192
  std::int32_t sector = 1;   // 1..12
  float r_cm = 23.2f;
  float phi_rad = 0.0f;
  float drift_time_ns = 0.0f;
  float pulse_height = 0.0f;
};

struct StwTrackElement {
  std::int32_t sector = 1;
  float r_ref_cm = 25.6f;
  float phi_rad = 0.0f;
  float dphi_dr = 0.0f;
  std::int32_t num_layers = 0;
};

[[nodiscard]] std::vector<StwHit> unpack_raw_blocklet(
    std::span<const std::int32_t> words);

[[nodiscard]] std::vector<StwTrackElement> reconstruct_elements(
    std::span<const StwHit> hits,
    float v_drift_cm_per_ns = 0.005f);

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const StwTrackElement> elements);

}  // namespace delphipp::reco::stw
