// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/calibration/global_alignment.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string_view>

#include "delphipp/reco/calibration/vd_alignment.hpp"

namespace delphipp::reco::calibration {

namespace {

std::string_view wire_db_path(WireDetector detector) noexcept {
  switch (detector) {
    case WireDetector::Id:
      return "GLOB$ID*.TIM0";
    case WireDetector::Od:
      return "GLOB$OD*.TIM0";
    case WireDetector::Fca:
      return "GLOB$FCA*.TIM0";
    case WireDetector::Fcb:
      return "GLOB$FCB*.TIM0";
  }
  return "GLOB$ID*.TIM0";
}

}  // namespace

DtransMatrix make_euler_transform(const Vec3& translation_cm,
                                  const Vec3& euler_deg) noexcept {
  constexpr double kDegToRad = std::numbers::pi_v<double> / 180.0;
  const double phi = static_cast<double>(euler_deg[0]) * kDegToRad;
  const double theta = static_cast<double>(euler_deg[1]) * kDegToRad;
  const double psi = static_cast<double>(euler_deg[2]) * kDegToRad;

  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  const double cthe = std::cos(theta);
  const double sthe = std::sin(theta);
  const double cpsi = std::cos(psi);
  const double spsi = std::sin(psi);

  // Matches glbali.car DXMATR Euler rotation convention.
  DtransMatrix rmatr{};
  rmatr[0] = translation_cm[0];
  rmatr[1] = translation_cm[1];
  rmatr[2] = translation_cm[2];

  rmatr[3] = static_cast<float>(cpsi * cphi - spsi * cthe * sphi);
  rmatr[4] = static_cast<float>(-spsi * cphi - cpsi * cthe * sphi);
  rmatr[5] = static_cast<float>(sthe * sphi);

  rmatr[6] = static_cast<float>(cpsi * sphi + spsi * cthe * cphi);
  rmatr[7] = static_cast<float>(-spsi * sphi + cpsi * cthe * cphi);
  rmatr[8] = static_cast<float>(-sthe * cphi);

  rmatr[9] = static_cast<float>(spsi * sthe);
  rmatr[10] = static_cast<float>(cpsi * sthe);
  rmatr[11] = static_cast<float>(cthe);
  return rmatr;
}

void GlobalDetectorAlignment::set_transform(std::int32_t detector_id,
                                            const Vec3& translation_cm,
                                            const Vec3& angles,
                                            bool euler_degrees) {
  transforms_[detector_id] = euler_degrees
                                 ? make_euler_transform(translation_cm, angles)
                                 : make_rigid_transform(translation_cm, angles);
}

void GlobalDetectorAlignment::set_matrix(std::int32_t detector_id,
                                         const DtransMatrix& rmatr) {
  transforms_[detector_id] = rmatr;
}

bool GlobalDetectorAlignment::has_transform(std::int32_t detector_id) const noexcept {
  return transforms_.find(detector_id) != transforms_.end();
}

std::optional<DtransMatrix> GlobalDetectorAlignment::matrix(
    std::int32_t detector_id) const noexcept {
  const auto it = transforms_.find(detector_id);
  if (it == transforms_.end()) {
    return std::nullopt;
  }
  return it->second;
}

Vec3 GlobalDetectorAlignment::apply_global_transform(const Vec3& point_cm,
                                                     std::int32_t detector_id,
                                                     bool inverse) const noexcept {
  const auto it = transforms_.find(detector_id);
  if (it == transforms_.end()) {
    return point_cm;
  }
  if (!inverse) {
    return delphipp::dtrans(point_cm, it->second, 0);
  }
  const DtransMatrix inv = invert_rigid_transform(it->second);
  return delphipp::dtrans(point_cm, inv, 0);
}

TpcCorrectionResult correct_tpc_coordinates_detailed(
    const Vec3& xyz_cm,
    const TpcCorrectionConfig& cfg,
    const db::ConditionDatabase* db,
    const MagneticField* bfield) noexcept {
  TpcCorrectionResult out;
  out.corrected_xyz_cm = xyz_cm;
  out.applied_abst0_bins = cfg.default_abst0_bins;

  const float x = xyz_cm[0];
  const float y = xyz_cm[1];
  const float z = xyz_cm[2];
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
    return out;
  }

  const double r = std::hypot(static_cast<double>(x), static_cast<double>(y));
  const double phi = (r > 1.0e-6) ? std::atan2(static_cast<double>(y),
                                               static_cast<double>(x))
                                  : 0.0;
  const std::string_view enp_path =
      (z >= 0.0f) ? "GLOB$TPC*.GLOB$ENP0" : "GLOB$TPC*.GLOB$ENP1";

  // 1. ConditionDatabase 94D ABST0 z-offset (tpct0_94c.asc: 1.147 -> 1.10208 => +-232 um).
  float v_drift = cfg.drift_velocity_cm_per_us;
  if (cfg.enable_abst0_correction) {
    if (db != nullptr) {
      if (const auto w14 = db->resolve_word(
              "CALB", enp_path, "SLOW", 14, cfg.query_time, cfg.entry_cutoff);
          w14.has_value() && std::isfinite(*w14)) {
        out.applied_abst0_bins = *w14;
      }
      if (cfg.use_db_drift_velocity) {
        if (const auto w1 = db->resolve_word(
                "CALB", enp_path, "SLOW", 1, cfg.query_time, cfg.entry_cutoff);
            w1.has_value() && std::isfinite(*w1) && *w1 > 0.0f) {
          v_drift = *w1;
        }
      }
    }
    const double delta_abs_z =
        static_cast<double>(out.applied_abst0_bins - cfg.reference_abst0_bins) *
        static_cast<double>(cfg.time_bin_us) * static_cast<double>(v_drift);
    out.delta_z_abst0_cm =
        static_cast<float>((z >= 0.0f) ? delta_abs_z : -delta_abs_z);
  }

  // 2. B-field homogeneity map correction & ExB radial shift.
  double bz_t = static_cast<double>(kSolenoidNominalT);
  double br_t = 0.0;
  double bphi_t = 0.0;
  if (bfield != nullptr) {
    const FieldVector fv = bfield->evaluate(xyz_cm);
    const double evaluated_bz = static_cast<double>(fv.bz_tesla());
    if (std::isfinite(evaluated_bz) && std::fabs(evaluated_bz) > 0.05) {
      bz_t = evaluated_bz;
      if (r > 1.0e-6) {
        const double bx_t = static_cast<double>(fv.bx) /
                            static_cast<double>(kFieldConvGeVPerCmT);
        const double by_t = static_cast<double>(fv.by) /
                            static_cast<double>(kFieldConvGeVPerCmT);
        const double c = std::cos(phi);
        const double s = std::sin(phi);
        br_t = bx_t * c + by_t * s;
        bphi_t = -bx_t * s + by_t * c;
      }
    }
  }

  const double wt = static_cast<double>(cfg.omega_tau) *
                    (bz_t / static_cast<double>(kSolenoidNominalT));
  const double denom = 1.0 + wt * wt;

  if (cfg.enable_exb_correction) {
    out.delta_r_exb_cm = static_cast<float>(
        static_cast<double>(cfg.exb_radial_shift_cm) * (wt / denom));
  }

  if (cfg.enable_bfield_map_correction && bfield != nullptr && std::fabs(bz_t) > 0.05) {
    const double l_drift = std::max(0.0, static_cast<double>(cfg.endplate_z_cm) -
                                             std::fabs(static_cast<double>(z)));
    const double sgn_z = (z >= 0.0f) ? 1.0 : -1.0;
    const double rel_br = br_t / bz_t;
    const double rel_bphi = bphi_t / bz_t;

    const double dr = (-(wt * wt / denom) * rel_br - (wt / denom) * rel_bphi) *
                      l_drift * sgn_z;
    const double drphi = (-(wt / denom) * rel_br + (wt * wt / denom) * rel_bphi) *
                         l_drift * sgn_z;

    out.delta_r_bfield_cm = static_cast<float>(dr);
    out.delta_rphi_bfield_cm = static_cast<float>(drphi);
  }

  const double r_new = std::max(
      0.0, r + static_cast<double>(out.delta_r_exb_cm) +
               static_cast<double>(out.delta_r_bfield_cm));
  const double dphi =
      (r > 1.0e-6) ? (static_cast<double>(out.delta_rphi_bfield_cm) / r) : 0.0;
  const double phi_new = phi + dphi;

  out.corrected_xyz_cm[0] = static_cast<float>(r_new * std::cos(phi_new));
  out.corrected_xyz_cm[1] = static_cast<float>(r_new * std::sin(phi_new));
  out.corrected_xyz_cm[2] = z + out.delta_z_abst0_cm;
  return out;
}

Vec3 correct_tpc_coordinates(const Vec3& xyz_cm,
                             const TpcCorrectionConfig& cfg,
                             const db::ConditionDatabase* db,
                             const MagneticField* bfield) noexcept {
  return correct_tpc_coordinates_detailed(xyz_cm, cfg, db, bfield).corrected_xyz_cm;
}

WireTimingConfig default_wire_timing_config(WireDetector detector,
                                            std::int32_t run_number) noexcept {
  const bool is_mc = (run_number <= 0);
  switch (detector) {
    case WireDetector::Id:
      return WireTimingConfig{
          .t0_ns = is_mc ? 0.0f : 12.5f,
          .t0_secondary_ns = 0.0f,
          .v_drift_cm_per_ns = 0.0040f,
          .v_prop_cm_per_ns = 20.0f,
          .z_scale = 1.0f,
          .z_offset_cm = is_mc ? 0.0f : 0.0060f,  // idbfix.car 93B/93C trigger z shift
          .rphi_offset_cm = 0.0f,
      };
    case WireDetector::Od:
      return WireTimingConfig{
          .t0_ns = is_mc ? 0.0f : 8.0f,
          .t0_secondary_ns = is_mc ? 0.0f : 8.0f,
          .v_drift_cm_per_ns = 0.0050f,
          .v_prop_cm_per_ns = 20.0f,
          .z_scale = is_mc ? 1.0033f : 1.0f,  // odbfix.car line 89 MC z scale
          .z_offset_cm = 0.0f,
          .rphi_offset_cm = 0.0f,
      };
    case WireDetector::Fca:
      return WireTimingConfig{
          .t0_ns = is_mc ? 0.0f : 15.0f,
          .t0_secondary_ns = 0.0f,
          .v_drift_cm_per_ns = 0.0050f,
          .v_prop_cm_per_ns = 20.0f,
          .z_scale = 1.0f,
          .z_offset_cm = 0.0f,
          .rphi_offset_cm = 0.0f,
      };
    case WireDetector::Fcb:
      return WireTimingConfig{
          .t0_ns = is_mc ? 0.0f : 153.0f,  // fcbt0.car TIMM nominal mean drift t0
          .t0_secondary_ns = 0.0f,
          .v_drift_cm_per_ns = 0.0050f,
          .v_prop_cm_per_ns = 20.0f,
          .z_scale = 1.0f,
          .z_offset_cm = 0.0f,
          .rphi_offset_cm = 0.0f,
      };
  }
  return WireTimingConfig{};
}

WireTimingConfig resolve_wire_timing_config(const db::ConditionDatabase& db,
                                            WireDetector detector,
                                            db::DbTimestamp query_time,
                                            std::int32_t run_number,
                                            db::DbTimestamp entry_cutoff) {
  WireTimingConfig cfg = default_wire_timing_config(detector, run_number);
  const auto words = db.resolve_field(
      "CALB", wire_db_path(detector), "CALI", query_time, entry_cutoff);
  if (const auto it = words.find(1); it != words.end() && std::isfinite(it->second)) {
    cfg.t0_ns = it->second;
    cfg.t0_secondary_ns = it->second;
  }
  if (const auto it = words.find(2);
      it != words.end() && std::isfinite(it->second) && it->second > 0.0f) {
    cfg.v_drift_cm_per_ns = it->second;
  }
  if (const auto it = words.find(3); it != words.end() && std::isfinite(it->second)) {
    cfg.z_offset_cm = it->second;
  }
  if (const auto it = words.find(4); it != words.end() && std::isfinite(it->second)) {
    cfg.rphi_offset_cm = it->second;
  }
  if (const auto it = words.find(5); it != words.end() && std::isfinite(it->second)) {
    cfg.t0_secondary_ns = it->second;
  }
  return cfg;
}

WireTimingResult apply_wire_t0_and_vdrift(const WireTimingInput& hit,
                                          const WireTimingConfig& cfg) noexcept {
  WireTimingResult res;
  if (!std::isfinite(hit.raw_time_ns) || !std::isfinite(cfg.t0_ns) ||
      !std::isfinite(cfg.v_drift_cm_per_ns) || cfg.v_drift_cm_per_ns <= 0.0f) {
    res.valid = false;
    return res;
  }

  const float t_corr = hit.raw_time_ns - cfg.t0_ns;
  res.corrected_time_ns = t_corr;
  res.valid = (t_corr >= 0.0f);
  res.drift_distance_cm = std::max(0.0f, t_corr * cfg.v_drift_cm_per_ns);
  const float sign = (hit.lr_sign >= 0) ? 1.0f : -1.0f;
  res.coordinate_cm =
      hit.wire_coord_cm + sign * res.drift_distance_cm + cfg.rphi_offset_cm;

  if (hit.has_dual_end_timing && hit.detector == WireDetector::Od &&
      std::isfinite(hit.secondary_time_ns)) {
    const float ta = hit.raw_time_ns - cfg.t0_ns;
    const float tc = hit.secondary_time_ns - cfg.t0_secondary_ns;
    const float raw_z = 0.5f * cfg.v_prop_cm_per_ns * (tc - ta);
    res.z_cm = raw_z * cfg.z_scale + cfg.z_offset_cm;
  } else {
    res.z_cm = hit.z_cm * cfg.z_scale + cfg.z_offset_cm;
  }
  return res;
}

}  // namespace delphipp::reco::calibration
