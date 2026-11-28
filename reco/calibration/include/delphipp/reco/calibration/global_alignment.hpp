// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI global detector alignment, TPC systematic
// coordinate/ABST0 fixers, and wire-chamber t0/vdrift calibrations
// (simana_v94d/src/car/glbali.car, tpcfix.car, idbfix.car, odbfix.car,
// odtiming.car, fcafix.car, fcbfix.car, fcbt0.car).

#pragma once

#include <cstdint>
#include <map>
#include <optional>

#include "delphipp/db/condition_db.hpp"
#include "delphipp/geometry.hpp"
#include "delphipp/magnetic_field.hpp"

namespace delphipp::reco::calibration {

// Standard DELPHI TANAGRA detector identifiers used in glbali.car.
enum class DetectorId : std::int32_t {
  Vd = 3,
  Id = 4,
  Tpc = 5,
  Brich = 7,
  Od = 8,
  Hpc = 9,
  Hac = 13,
  Mub = 14,
  Mus = 17,
  Fca = 20,
  Fcb = 25,
  Emf = 26,
  Frich = 27,
  Muf = 30
};

// Build DDAPP transformation matrix from translation [cm] and Euler angles [deg] (glbali.car DXMATR).
[[nodiscard]] DtransMatrix make_euler_transform(
    const Vec3& translation_cm,
    const Vec3& euler_deg) noexcept;

// Per-detector global alignment store matching glbali.car ALMINI / TRKROT.
class GlobalDetectorAlignment {
 public:
  void set_transform(std::int32_t detector_id,
                     const Vec3& translation_cm,
                     const Vec3& angles,
                     bool euler_degrees = false);

  void set_transform(DetectorId detector_id,
                     const Vec3& translation_cm,
                     const Vec3& angles,
                     bool euler_degrees = false) {
    set_transform(static_cast<std::int32_t>(detector_id),
                  translation_cm,
                  angles,
                  euler_degrees);
  }

  void set_matrix(std::int32_t detector_id, const DtransMatrix& rmatr);

  [[nodiscard]] bool has_transform(std::int32_t detector_id) const noexcept;
  [[nodiscard]] std::optional<DtransMatrix> matrix(
      std::int32_t detector_id) const noexcept;

  // Apply forward (local->DELPHI) or inverse (DELPHI->local) 6-DOF transform via delphipp::dtrans.
  [[nodiscard]] Vec3 apply_global_transform(const Vec3& point_cm,
                                            std::int32_t detector_id,
                                            bool inverse = false) const noexcept;

  [[nodiscard]] Vec3 apply_global_transform(const Vec3& point_cm,
                                            DetectorId detector_id,
                                            bool inverse = false) const noexcept {
    return apply_global_transform(
        point_cm, static_cast<std::int32_t>(detector_id), inverse);
  }

 private:
  std::map<std::int32_t, DtransMatrix> transforms_;
};

[[nodiscard]] inline Vec3 apply_global_transform(
    const Vec3& point_cm,
    std::int32_t detector_id,
    const GlobalDetectorAlignment& alignment,
    bool inverse = false) noexcept {
  return alignment.apply_global_transform(point_cm, detector_id, inverse);
}

// DELPHI TPC constants (tpcana.car / ddb/tpct0_94c.asc).
inline constexpr float kTpcEndplateZCm = 134.0f;
inline constexpr float kTpcTimeBinUs = 26.0f / 352.209f;  // TPIX in tpcana.car
inline constexpr float kTpcNominalDriftVelCmPerUs = 7.00f;
inline constexpr float kTpcAbst0Sim = 1.147f;
inline constexpr float kTpcAbst0_94D = 1.10208f;

struct TpcCorrectionConfig {
  db::DbTimestamp query_time{940526, 120000};
  db::DbTimestamp entry_cutoff{0, 0};  // Default DELANA year-2000 cutoff
  float reference_abst0_bins = kTpcAbst0Sim;
  float default_abst0_bins = kTpcAbst0Sim;
  float drift_velocity_cm_per_us = kTpcNominalDriftVelCmPerUs;
  bool use_db_drift_velocity = false;
  float time_bin_us = kTpcTimeBinUs;
  float endplate_z_cm = kTpcEndplateZCm;
  float omega_tau = 1.8f;              // Electron Hall parameter mu * Bz at 1.2 T
  float exb_radial_shift_cm = 0.015f;  // Wire-plane ExB radial shift scale [cm]
  bool enable_exb_correction = true;
  bool enable_bfield_map_correction = true;
  bool enable_abst0_correction = true;
};

struct TpcCorrectionResult {
  Vec3 corrected_xyz_cm{};
  float delta_r_exb_cm = 0.0f;
  float delta_r_bfield_cm = 0.0f;
  float delta_rphi_bfield_cm = 0.0f;
  float delta_z_abst0_cm = 0.0f;
  float applied_abst0_bins = kTpcAbst0Sim;
};

// Apply TPC ExB radial shift, B-field homogeneity map correction, and ConditionDatabase
// 94D ABST0 z-offset (tpcfix.car + ddb/tpct0_94c.asc).
[[nodiscard]] TpcCorrectionResult correct_tpc_coordinates_detailed(
    const Vec3& xyz_cm,
    const TpcCorrectionConfig& cfg = {},
    const db::ConditionDatabase* db = nullptr,
    const MagneticField* bfield = nullptr) noexcept;

[[nodiscard]] Vec3 correct_tpc_coordinates(
    const Vec3& xyz_cm,
    const TpcCorrectionConfig& cfg = {},
    const db::ConditionDatabase* db = nullptr,
    const MagneticField* bfield = nullptr) noexcept;

// Wire chambers covered by t0/vdrift fixers (idbfix, odbfix, odtiming, fcafix, fcbfix, fcbt0).
enum class WireDetector : std::uint8_t {
  Id,
  Od,
  Fca,
  Fcb
};

struct WireTimingInput {
  WireDetector detector = WireDetector::Id;
  std::int32_t sector_or_module = 1;
  std::int32_t endcap_or_side = 1;  // 1: +z / C-side, 2: -z / A-side
  float raw_time_ns = 0.0f;
  float secondary_time_ns = 0.0f;   // Dual-end timing second end (OD C-end)
  float wire_coord_cm = 0.0f;
  float z_cm = 0.0f;
  std::int32_t lr_sign = 1;
  bool has_dual_end_timing = false;
};

struct WireTimingConfig {
  float t0_ns = 0.0f;
  float t0_secondary_ns = 0.0f;
  float v_drift_cm_per_ns = 0.0040f;
  float v_prop_cm_per_ns = 20.0f;
  float z_scale = 1.0f;
  float z_offset_cm = 0.0f;
  float rphi_offset_cm = 0.0f;
};

struct WireTimingResult {
  float corrected_time_ns = 0.0f;
  float drift_distance_cm = 0.0f;
  float coordinate_cm = 0.0f;
  float z_cm = 0.0f;
  bool valid = true;
};

// Nominal run-dependent wire timing defaults (ID, OD, FCA, FCB).
[[nodiscard]] WireTimingConfig default_wire_timing_config(
    WireDetector detector,
    std::int32_t run_number = 47014) noexcept;

// Resolve wire timing constants from ConditionDatabase (table "CALB", field "CALI"),
// falling back to default_wire_timing_config when absent.
[[nodiscard]] WireTimingConfig resolve_wire_timing_config(
    const db::ConditionDatabase& db,
    WireDetector detector,
    db::DbTimestamp query_time,
    std::int32_t run_number = 47014,
    db::DbTimestamp entry_cutoff = db::DbTimestamp{0, 0});

// Apply t0 subtraction, drift-velocity conversion, dual-end z timing (OD), and z/rphi offsets.
[[nodiscard]] WireTimingResult apply_wire_t0_and_vdrift(
    const WireTimingInput& hit,
    const WireTimingConfig& cfg) noexcept;

}  // namespace delphipp::reco::calibration
