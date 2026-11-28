// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Silicon Vertex Detector (VD) offline alignment,
// fill-by-fill beamspot fitting, and pedestal/gain calibration
// (simana_v94d/src/car/vdalign.car, vdbeam.car, vdcalyd.car).

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "delphipp/geometry.hpp"

namespace delphipp::reco::calibration {

// 6-DOF rigid survey transform for a VD ladder/plaquette (vdalign.car / glbali.car DCMATR).
struct VdLadderSurvey {
  std::int32_t layer = 1;   // 1=Closer, 2=Inner, 3=Outer
  std::int32_t module = 1;  // Ladder/module ID (e.g., 1..24 or +-layer*100+sector)
  Vec3 translation_cm{};    // Local-to-DELPHI translation (dx, dy, dz) [cm]
  Vec3 rotation_rad{};      // Cartesian rotation angles (rx, ry, rz) [rad]
  float humidity_bow_cm = 0.0f;  // Optional radial bow sagitta [cm] (vdalign.car EVHUBE)

  [[nodiscard]] DtransMatrix to_dtrans_matrix() const noexcept;
  [[nodiscard]] DtransMatrix to_inverse_dtrans_matrix() const noexcept;
};

// Build a 12-word DDAPP DTRANS matrix from translation [cm] and Cartesian angles [rad] (DCMATR).
[[nodiscard]] DtransMatrix make_rigid_transform(
    const Vec3& translation_cm,
    const Vec3& rotation_rad) noexcept;

// Invert a 12-word DDAPP DTRANS rigid transform (ALINVE).
[[nodiscard]] DtransMatrix invert_rigid_transform(
    std::span<const float, 12> rmatr) noexcept;

// Transform a VD point between local ladder coordinates and the DELPHI global frame
// using 6-DOF survey rigid transforms via delphipp::dtrans.
// If no matching (layer, module) survey is found, returns `point_cm` unchanged.
[[nodiscard]] Vec3 align_vd_point(
    const Vec3& point_cm,
    std::int32_t layer,
    std::int32_t module,
    std::span<const VdLadderSurvey> surveys,
    bool to_global = true) noexcept;

// Convenience overload for a single survey entry.
[[nodiscard]] Vec3 align_vd_point(
    const Vec3& point_cm,
    std::int32_t layer,
    std::int32_t module,
    const VdLadderSurvey& survey,
    bool to_global = true) noexcept;

// Configuration for iterative truncated beamspot fit (vdbeam.car).
struct BeamspotFitConfig {
  float sigma_cut = 3.0f;
  std::size_t max_iterations = 12;
  Vec3 min_scale_cm = {0.0010f, 0.0010f, 0.0100f};  // 10 um transverse, 100 um z floor
  double convergence_tol_cm = 1e-6;
};

// Result of a fill-by-fill iterative truncated primary-vertex beamspot fit.
struct BeamspotFitResult {
  std::int32_t fill_number = 0;
  Vec3 position_cm{};
  Vec3 uncertainty_cm{};   // Standard error on beamspot centroid [cm]
  Vec3 beam_size_cm{};     // Fitted RMS beam envelope (sigma_x, sigma_y, sigma_z) [cm]
  std::array<float, 6> covariance_cm2{};  // Symmetric 3x3: [cxx, cxy, cxz, cyy, cyz, czz]
  std::size_t input_count = 0;
  std::size_t accepted_count = 0;
  std::size_t iterations = 0;
  bool converged = false;
  std::int32_t quality_flag = 2;  // vdbeam.car IERR: 0=good, 1=low-stat, 2=failed
};

// Robust iterative truncated mean & covariance fit rejecting displaced secondary vertices.
[[nodiscard]] BeamspotFitResult fit_fill_beamspot(
    std::span<const Vec3> primary_vertices,
    std::int32_t fill_number,
    const BeamspotFitConfig& cfg = {});

// Channel quality status codes matching vdcalyd.car REDATA/REACAL.
enum class VdChannelStatus : std::uint8_t {
  Good = 0,
  Dead = 1,
  Grey = 5,
  HotNoisy = 6,
  BadGain = 7,
  PunchThrough = 8
};

// Per-strip calibration record for pedestal subtraction, gain equalization, and masking.
struct VdChannelCalibration {
  std::int32_t layer = 1;
  std::int32_t module = 1;
  std::int32_t strip = 0;
  float pedestal = 0.0f;
  float relative_gain = 1.0f;  // Calibrated signal = (raw_adc - pedestal) / relative_gain
  float noise = 1.0f;
  VdChannelStatus status = VdChannelStatus::Good;
};

struct VdRawStripInput {
  std::int32_t layer = 1;
  std::int32_t module = 1;
  std::int32_t strip = 0;
  float raw_adc = 0.0f;
  bool is_z_side = false;
};

struct VdCalibratedStrip {
  std::int32_t layer = 1;
  std::int32_t module = 1;
  std::int32_t strip = 0;
  float pulse_height = 0.0f;
  float noise = 1.0f;
  float sn_ratio = 0.0f;
  bool is_z_side = false;
};

// Mask dead/hot/punch-through channels and apply pedestal subtraction + gain equalization.
// Returns std::nullopt if the channel is masked or has invalid calibration.
[[nodiscard]] std::optional<float> apply_vd_pedestal_gain(
    float raw_adc,
    const VdChannelCalibration& calib) noexcept;

// Batch calibration for raw VD strips using a channel calibration lookup table.
[[nodiscard]] std::vector<VdCalibratedStrip> apply_vd_pedestal_gain(
    std::span<const VdRawStripInput> raw_strips,
    std::span<const VdChannelCalibration> calib_table);

}  // namespace delphipp::reco::calibration
