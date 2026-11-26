// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Small Angle Tile Calorimeter (STIC / STC, 1994+
// primary luminosity monitor, simana_v94d/src/car/stcana.car).

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::stc {

inline constexpr std::int32_t kDetectorId = 32;  // STIC detector code
inline constexpr std::int32_t kNumRings = 10;    // 10 radial rings
inline constexpr std::int32_t kNumSectors = 16;  // 16 azimuthal sectors per arm
inline constexpr float kNominalZCm = 220.0f;     // Front face |z| [cm]

struct StcTowerHit {
  std::int32_t arm = 1;       // 1 (+z) or 2 (-z)
  std::int32_t ring = 1;      // 1..10
  std::int32_t sector = 1;    // 1..16 (cyclic)
  float energy_gev = 0.0f;
  float r_cm = 15.0f;
  float phi_rad = 0.0f;
  float z_cm = kNominalZCm;
};

struct StcCluster {
  std::int32_t arm = 1;
  float energy_gev = 0.0f;
  float r_cm = 0.0f;
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
  std::int32_t num_towers = 0;
};

[[nodiscard]] std::vector<StcTowerHit> unpack_raw_blocklet(
    std::span<const std::int32_t> words);

[[nodiscard]] std::vector<StcCluster> find_clusters(
    std::span<const StcTowerHit> towers,
    float min_cluster_energy_gev = 0.5f);

// Basic back-to-back Bhabha luminosity event check between Arm 1 and Arm 2.
[[nodiscard]] bool is_bhabha_candidate(
    std::span<const StcCluster> clusters,
    float min_arm_energy_gev = 25.0f,
    float max_acoplanarity_rad = 0.35f) noexcept;

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const StcCluster> clusters);

}  // namespace delphipp::reco::stc
