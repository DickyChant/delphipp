// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Surround Muon Chambers (MUS) reconstruction
// (simana_v94d/src/car/musana.car: EDAMSH, EDCLSR, EDCBCL, EDBLTE).

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/geometry.hpp"

namespace delphipp::reco::mus {

inline constexpr std::int32_t kDetectorId = 17;  // MUS detector code
inline constexpr std::size_t kNumDetectors = 16;
inline constexpr std::size_t kNumLayers = 4;     // 1:strip, 2:tube, 3:tube, 4:strip
inline constexpr std::size_t kMaxHitsPerLayer = 50;
inline constexpr std::size_t kMaxClustersPerLayer = 6;
inline constexpr std::size_t kMaxHitsPerCluster = 15;
inline constexpr std::size_t kMaxSpacePoints = 20;
inline constexpr std::size_t kMaxSolutionsPerPoint = 5;

struct MusConfig {
  std::int32_t strip_pedestal = 30;  // IPEDS (musana.car 19.11.96 AO)
  std::int32_t anode_pedestal = 16;  // IPEDA
  float sigma_strip_cm = 1.0f;       // SIGSTR
  float sigma_tube_cm = 2.5f;        // SIGTUB
  float beta_bottom_rad = 0.12f;     // 6.9 deg
  float beta_top_side_rad = 0.84f;   // 48 deg
};

// Geometry calibration constants from /EDGDDB/ (XYZTS1, DXYTS1, TMPLLD)
struct MusGeometry {
  // [det 0..15][layer 0..3] -> origin (x, y, z) of element 1 in local frame
  std::array<std::array<Vec3, kNumLayers>, kNumDetectors> origin_ts1{};
  // [det 0..15][layer 0..3] -> pitch (dx, dy)
  std::array<std::array<std::array<float, 2>, kNumLayers>, kNumDetectors> pitch_ts1{};
  // [det 0..15][plane 0..1] -> local-to-DELPHI DTRANS 12-word affine matrix
  std::array<std::array<DtransMatrix, 2>, kNumDetectors> local_to_delphi{};
};

struct MusRawHit {
  std::int32_t det = 0;      // 1..16
  std::int32_t layer = 0;    // 1..4 (1,4 = strips; 2,3 = tubes)
  std::int32_t element = 0;  // 1-based tube or strip number (ITSNB)
  float amplitude = 0.0f;    // calibrated amplitude (AMPL)
};

struct MusCluster {
  std::int32_t det = 0;
  std::int32_t layer = 0;
  float center_of_mass = 0.0f;  // amplitude-weighted element number (CMCLU)
  float total_charge = 0.0f;    // sum of amplitudes (CHCLU)
  std::vector<std::size_t> raw_hit_indices;
};

struct MusSpacePointSolution {
  Vec3 local_xyz{};
  std::size_t tube_cluster_idx = 0;
  std::size_t strip_cluster_idx = 0;
};

struct MusSpacePoint {
  std::int32_t det = 0;        // 1..16
  std::int32_t plane = 0;      // 1 or 2
  std::int32_t submod_id = 0;  // det * 10 + plane
  std::vector<MusSpacePointSolution> solutions;
};

// 20-word MUS TE bank payload in cylindrical coordinates (R, R*phi, z)
struct MusTeBank {
  std::int32_t detector_id = kDetectorId;
  std::int32_t submod_id = 0;
  std::int32_t descriptor = 125;  // (R, Rphi, z)
  float r = 0.0f;
  float r_phi = 0.0f;
  float z = 0.0f;
  float cov_rphi = 0.0f;     // (SIGSTR * cos(alpha))^2
  float cov_rphi_z = 0.0f;   // 0
  float cov_z = 0.0f;        // (SIGTUB / cos(beta))^2
};

// Port of SUBROUTINE EDAMSH: returns -1.0f if below threshold, else float(adc).
[[nodiscard]] float shape_amplitude(std::int32_t iadr,
                                    std::int32_t iadc,
                                    const MusConfig& cfg = {}) noexcept;

// Port of SUBROUTINE EDCLSR: cluster contiguous element hits per (det, layer).
[[nodiscard]] std::vector<MusCluster> find_clusters(std::span<const MusRawHit> raw_hits);

// Port of SUBROUTINE EDCBCL: combine tube clusters (layers 2,3) with complementary
// strip clusters (layers 1,4) into local space points (preserving ambiguities).
[[nodiscard]] std::vector<MusSpacePoint> combine_clusters(
    std::span<const MusCluster> clusters,
    const MusGeometry& geom);

// Port of SUBROUTINE EDBLTE: transform space points to DELPHI cylindrical coords
// and compute TE covariance matrices.
[[nodiscard]] std::vector<MusTeBank> build_te(
    std::span<const MusSpacePoint> space_points,
    const MusGeometry& geom,
    const MusConfig& cfg = {});

}  // namespace delphipp::reco::mus
