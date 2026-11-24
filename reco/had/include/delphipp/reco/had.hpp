// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Hadron Calorimeter (HAC / HCAL) reconstruction
// (simana_v94d/src/car/hadana.car: ECALIB tower calibration with pressure,
// isobutane, and temperature corrections, plus 3D tower shower clustering).

#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::had {

inline constexpr std::int32_t kDetectorId = 13;  // HAC detector code in TANAGRA (barrel/endcap)
inline constexpr std::int32_t kNumPhiBins = 96;  // Cyclic phi segmentation

// Mirrors ECGASD gas/environment parameters used in ECALIB (hadana.car lines 2049-2061)
struct HadEnvironment {
  std::int32_t run_number = 47014;  // Real data if > 10581
  float pressure_mbar = 975.0f;     // ECGASD(2): active in (950, 1000)
  float isobutane_pct = 30.0f;      // ECGASD(7): active in (29, 31)
  float temperature_c = 20.0f;      // ECGASD(11): active in (16, 23)
  float global_cal_factor = 1.0f;   // CLBFCR
};

struct HadTowerHit {
  std::int32_t itheta = 0;  // Polar tower index
  std::int32_t iphi = 0;    // Azimuthal tower index (0..95, cyclic)
  std::int32_t ilayer = 0;  // Depth layer (1..4)
  float energy_gev = 0.0f;
  float r_cm = 0.0f;
  float theta_rad = 0.0f;
  float phi_rad = 0.0f;
};

struct HadCluster {
  float energy_gev = 0.0f;
  float r_cm = 0.0f;
  float theta_rad = 0.0f;
  float phi_rad = 0.0f;
  std::int32_t num_towers = 0;
  std::int32_t max_layer = 0;
};

// Port of SUBROUTINE ECALIB (hadana.car lines 2043-2069).
// Returns calibrated tower energy in GeV, or std::nullopt if below threshold.
[[nodiscard]] std::optional<float> calibrate_tower(
    std::int32_t signal,
    std::int32_t threshold,
    std::int32_t tower_coeff_mev,
    const HadEnvironment& env = {}) noexcept;

// 3D (itheta, cyclic iphi, ilayer) connected-component tower clustering.
[[nodiscard]] std::vector<HadCluster> find_clusters(
    std::span<const HadTowerHit> towers,
    float min_cluster_energy_gev = 0.3f);

// Write reconstructed HAC clusters as TER banks into TanagraStore.
std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const HadCluster> clusters);

}  // namespace delphipp::reco::had
