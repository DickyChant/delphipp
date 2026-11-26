// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI High-Density Projection Chamber (HPC / Barrel ECAL)
// reconstruction (simana_v94d/src/car/hpcana.car: drift-to-z conversion, 3D
// shower clustering across 9 longitudinal layers, and TANAGRA TER export).

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::hpc {

inline constexpr std::int32_t kDetectorId = 9;     // HPC detector code
inline constexpr std::size_t kNumLayers = 9;       // 9 longitudinal sampling layers
inline constexpr float kInnerRadiusCm = 208.0f;    // Nominal HPC inner radius [cm]
inline constexpr float kOuterRadiusCm = 260.0f;    // Nominal HPC outer radius [cm]

struct HpcConfig {
  float drift_velocity_cm_per_us = 5.4f;
  float min_charge_gev = 0.02f;
  float min_shower_energy_gev = 0.15f;
  float cluster_radius_rphi_cm = 8.0f;
  float cluster_radius_z_cm = 8.0f;
};

struct HpcChargeHit {
  std::int32_t module = 1;       // 1..144 (24 phi * 6 z rings)
  std::int32_t layer = 1;        // 1..9
  float r_cm = 225.0f;
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
  float energy_gev = 0.0f;
};

struct HpcShower {
  std::int32_t module = 1;
  float energy_gev = 0.0f;
  float r_cm = 0.0f;
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
  std::array<float, kNumLayers> layer_energy_gev{};
  std::int32_t num_hits = 0;
};

[[nodiscard]] std::vector<HpcChargeHit> unpack_raw_blocklet(
    std::span<const std::int32_t> words,
    const HpcConfig& cfg = {});

[[nodiscard]] std::vector<HpcShower> reconstruct_showers(
    std::span<const HpcChargeHit> hits,
    const HpcConfig& cfg = {});

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const HpcShower> showers);

}  // namespace delphipp::reco::hpc
