// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Forward Electromagnetic Calorimeter (FEMC / EMF)
// reconstruction (simana_v94d/src/car/emfana.car: calibration, EEFIND 2D block
// shower clustering, centroid & dispersion, TANAGRA TER export).

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::emf {

inline constexpr std::int32_t kDetectorId = 26;     // EMF detector code
inline constexpr float kNominalZCm = 284.0f;        // Nominal FEMC front face |z| [cm]
inline constexpr float kBlockPitchCm = 5.6f;        // Lead-glass block transverse pitch [cm]

struct EmfConfig {
  float min_block_energy_gev = 0.05f;   // Minimum block energy for clustering
  float min_cluster_energy_gev = 0.20f; // Minimum shower energy to keep
  float ph_per_gev_mc = 1340.0f;        // PHGEV in EECNT (emfana.car line 774)
  float max_cluster_energy_gev = 150.0f;// TABEMX
};

struct EmfBlockHit {
  std::int32_t endcap = 1;  // 1 (+z, side C/A) or 2 (-z)
  std::int32_t ix = 0;      // Grid index X
  std::int32_t iy = 0;      // Grid index Y
  float energy_gev = 0.0f;  // Calibrated block energy
  float x_cm = 0.0f;        // Block center X in DELPHI frame
  float y_cm = 0.0f;        // Block center Y in DELPHI frame
  float z_cm = kNominalZCm; // Block Z in DELPHI frame
};

struct EmfCluster {
  std::int32_t endcap = 1;
  float energy_gev = 0.0f;
  float x_cm = 0.0f;
  float y_cm = 0.0f;
  float z_cm = 0.0f;
  float sigma_x_cm = 0.0f;
  float sigma_y_cm = 0.0f;
  std::int32_t num_blocks = 0;
};

// 8-connected 2D block clustering per endcap (matching EEFIND neighborhood logic).
[[nodiscard]] std::vector<EmfCluster> find_clusters(
    std::span<const EmfBlockHit> blocks,
    const EmfConfig& cfg = {});

// Write reconstructed EMF clusters as TER banks into TanagraStore.
std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const EmfCluster> clusters);

}  // namespace delphipp::reco::emf
