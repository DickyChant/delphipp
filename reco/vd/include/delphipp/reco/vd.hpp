// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Silicon Vertex Detector (VD) reconstruction
// (simana_v94d/src/car/vdana.car: 3-layer Closer/Inner/Outer silicon strip
// clustering with S/N threshold and charge-weighted centroid).

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::vd {

inline constexpr std::int32_t kDetectorId = 3;  // VD detector code
inline constexpr std::array<float, 3> kLayerRadiiCm = {6.3f, 9.0f, 10.9f};  // Closer, Inner, Outer

struct VdConfig {
  float seed_sn_cut = 4.0f;      // Minimum seed strip S/N
  float neighbor_sn_cut = 2.0f;  // Minimum neighbor strip S/N
  float strip_pitch_cm = 0.0050f;// 50 um readout pitch
};

struct VdStripHit {
  std::int32_t layer = 1;    // 1=Closer, 2=Inner, 3=Outer
  std::int32_t module = 1;   // Azimuthal ladder 1..24
  std::int32_t strip = 0;    // Strip index
  float pulse_height = 0.0f;
  float noise = 1.0f;
  bool is_z_side = false;    // false: R-phi side, true: z-side (1994+ double-sided)
};

struct VdCluster {
  std::int32_t layer = 1;
  std::int32_t module = 1;
  bool is_z_side = false;
  float centroid_strip = 0.0f;
  float local_coord_cm = 0.0f;
  float total_charge = 0.0f;
  float sn_ratio = 0.0f;
  std::int32_t width = 0;
};

[[nodiscard]] std::vector<VdCluster> find_clusters(
    std::span<const VdStripHit> strips,
    const VdConfig& cfg = {});

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const VdCluster> clusters);

}  // namespace delphipp::reco::vd
