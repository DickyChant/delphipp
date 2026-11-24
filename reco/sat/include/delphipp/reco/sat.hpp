// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Small Angle Tagger & VSAT Luminosity monitors
// (simana_v94d/src/car/satana.car, lumana.car, stwana.car).

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::sat {

inline constexpr std::int32_t kDetectorIdSat = 21;   // SAT
inline constexpr std::int32_t kDetectorIdVsat = 22;  // VSAT / LUM

struct SatElementHit {
  std::int32_t arm = 1;       // 1 (+z) or 2 (-z)
  std::int32_t ring = 1;      // 1..8
  std::int32_t sector = 1;    // 1..24
  float energy_gev = 0.0f;
  float r_cm = 20.0f;
  float phi_rad = 0.0f;
  float z_cm = 232.0f;
};

struct SatCluster {
  std::int32_t arm = 1;
  float energy_gev = 0.0f;
  float r_cm = 0.0f;
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
  std::int32_t num_hits = 0;
};

[[nodiscard]] std::vector<SatCluster> reconstruct_clusters(
    std::span<const SatElementHit> hits,
    float min_energy_gev = 1.0f);

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const SatCluster> clusters,
    std::int32_t detector_id = kDetectorIdSat);

}  // namespace delphipp::reco::sat
