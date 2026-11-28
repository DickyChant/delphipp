// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/calibration/vd_alignment.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace delphipp::reco::calibration {

namespace {

double compute_median(std::vector<double>& vals) {
  const std::size_t n = vals.size();
  const std::size_t mid = n / 2;
  std::nth_element(vals.begin(), vals.begin() + static_cast<std::ptrdiff_t>(mid), vals.end());
  if ((n % 2U) == 1U) {
    return vals[mid];
  }
  const double hi = vals[mid];
  const auto max_lo = std::max_element(
      vals.begin(), vals.begin() + static_cast<std::ptrdiff_t>(mid));
  return 0.5 * (*max_lo + hi);
}

bool is_finite_vec3(const Vec3& v) noexcept {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

bool is_masked_channel(VdChannelStatus status) noexcept {
  return status == VdChannelStatus::Dead ||
         status == VdChannelStatus::HotNoisy ||
         status == VdChannelStatus::PunchThrough;
}

}  // namespace

DtransMatrix make_rigid_transform(const Vec3& translation_cm,
                                  const Vec3& rotation_rad) noexcept {
  const double cx = std::cos(static_cast<double>(rotation_rad[0]));
  const double sx = std::sin(static_cast<double>(rotation_rad[0]));
  const double cy = std::cos(static_cast<double>(rotation_rad[1]));
  const double sy = std::sin(static_cast<double>(rotation_rad[1]));
  const double cz = std::cos(static_cast<double>(rotation_rad[2]));
  const double sz = std::sin(static_cast<double>(rotation_rad[2]));

  // Matches glbali.car DCMATR row-major 3x3 rotation convention.
  DtransMatrix rmatr{};
  rmatr[0] = translation_cm[0];
  rmatr[1] = translation_cm[1];
  rmatr[2] = translation_cm[2];

  rmatr[3] = static_cast<float>(cy * cz);
  rmatr[4] = static_cast<float>(-cx * sz + sx * sy * cz);
  rmatr[5] = static_cast<float>(sx * sz + cx * sy * cz);

  rmatr[6] = static_cast<float>(cy * sz);
  rmatr[7] = static_cast<float>(cx * cz + sx * sy * sz);
  rmatr[8] = static_cast<float>(-sx * cz + cx * sy * sz);

  rmatr[9] = static_cast<float>(-sy);
  rmatr[10] = static_cast<float>(sx * cy);
  rmatr[11] = static_cast<float>(cx * cy);
  return rmatr;
}

DtransMatrix invert_rigid_transform(std::span<const float, 12> rmatr) noexcept {
  // Matches glbali.car ALINVE: R_inv = R^T, t_inv = -R^T * t.
  DtransMatrix inv{};
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      inv[3 + 3 * i + j] = rmatr[3 + 3 * j + i];
    }
  }
  const Vec3 t{rmatr[0], rmatr[1], rmatr[2]};
  const Vec3 rotated_t = delphipp::dtrans(t, inv, 1);
  inv[0] = -rotated_t[0];
  inv[1] = -rotated_t[1];
  inv[2] = -rotated_t[2];
  return inv;
}

DtransMatrix VdLadderSurvey::to_dtrans_matrix() const noexcept {
  return make_rigid_transform(translation_cm, rotation_rad);
}

DtransMatrix VdLadderSurvey::to_inverse_dtrans_matrix() const noexcept {
  const auto fwd = to_dtrans_matrix();
  return invert_rigid_transform(fwd);
}

Vec3 align_vd_point(const Vec3& point_cm,
                    std::int32_t layer,
                    std::int32_t module,
                    std::span<const VdLadderSurvey> surveys,
                    bool to_global) noexcept {
  for (const auto& s : surveys) {
    if (s.layer == layer && s.module == module) {
      return align_vd_point(point_cm, layer, module, s, to_global);
    }
  }
  return point_cm;
}

Vec3 align_vd_point(const Vec3& point_cm,
                    std::int32_t layer,
                    std::int32_t module,
                    const VdLadderSurvey& survey,
                    bool to_global) noexcept {
  if (survey.layer != layer || survey.module != module) {
    return point_cm;
  }
  Vec3 adjusted = point_cm;
  if (to_global) {
    // Apply radial humidity sagitta along local y (outward normal) before global transform.
    adjusted[1] += survey.humidity_bow_cm;
    const DtransMatrix mat = survey.to_dtrans_matrix();
    return delphipp::dtrans(adjusted, mat, 0);
  }
  const DtransMatrix inv = survey.to_inverse_dtrans_matrix();
  Vec3 local = delphipp::dtrans(adjusted, inv, 0);
  local[1] -= survey.humidity_bow_cm;
  return local;
}

BeamspotFitResult fit_fill_beamspot(std::span<const Vec3> primary_vertices,
                                    std::int32_t fill_number,
                                    const BeamspotFitConfig& cfg) {
  BeamspotFitResult result;
  result.fill_number = fill_number;
  result.input_count = primary_vertices.size();

  std::vector<Vec3> valid;
  valid.reserve(primary_vertices.size());
  for (const auto& v : primary_vertices) {
    if (is_finite_vec3(v)) {
      valid.push_back(v);
    }
  }
  if (valid.empty()) {
    return result;
  }

  // Robust initial center and spread via component-wise median & MAD.
  std::array<double, 3> center{};
  std::array<double, 3> scale{};
  std::vector<double> scratch;
  scratch.reserve(valid.size());

  for (std::size_t k = 0; k < 3; ++k) {
    scratch.clear();
    for (const auto& v : valid) {
      scratch.push_back(static_cast<double>(v[k]));
    }
    center[k] = compute_median(scratch);

    scratch.clear();
    for (const auto& v : valid) {
      scratch.push_back(std::fabs(static_cast<double>(v[k]) - center[k]));
    }
    const double mad = 1.4826 * compute_median(scratch);
    scale[k] = std::max(mad, static_cast<double>(cfg.min_scale_cm[k]));
  }

  const double cut = std::max(1.0, static_cast<double>(cfg.sigma_cut));
  const double chi2_cut = cut * cut * 3.0;  // 3D normalized distance threshold

  std::vector<std::uint8_t> active(valid.size(), 1U);
  std::array<double, 6> cov{};
  std::size_t kept = 0;

  for (std::size_t iter = 0; iter < std::max<std::size_t>(1U, cfg.max_iterations); ++iter) {
    result.iterations = iter + 1U;
    std::vector<std::uint8_t> next_active(valid.size(), 0U);
    std::size_t next_kept = 0;

    for (std::size_t i = 0; i < valid.size(); ++i) {
      double norm2 = 0.0;
      bool axis_ok = true;
      for (std::size_t k = 0; k < 3; ++k) {
        const double d = (static_cast<double>(valid[i][k]) - center[k]) / scale[k];
        if (std::fabs(d) > cut) {
          axis_ok = false;
        }
        norm2 += d * d;
      }
      if (axis_ok && norm2 <= chi2_cut) {
        next_active[i] = 1U;
        ++next_kept;
      }
    }

    if (next_kept == 0U) {
      break;
    }

    std::array<double, 3> new_center{};
    for (std::size_t i = 0; i < valid.size(); ++i) {
      if (next_active[i] != 0U) {
        for (std::size_t k = 0; k < 3; ++k) {
          new_center[k] += static_cast<double>(valid[i][k]);
        }
      }
    }
    const double inv_n = 1.0 / static_cast<double>(next_kept);
    for (std::size_t k = 0; k < 3; ++k) {
      new_center[k] *= inv_n;
    }

    std::array<double, 6> new_cov{};
    if (next_kept > 1U) {
      const double denom = static_cast<double>(next_kept - 1U);
      for (std::size_t i = 0; i < valid.size(); ++i) {
        if (next_active[i] == 0U) {
          continue;
        }
        const double dx = static_cast<double>(valid[i][0]) - new_center[0];
        const double dy = static_cast<double>(valid[i][1]) - new_center[1];
        const double dz = static_cast<double>(valid[i][2]) - new_center[2];
        new_cov[0] += dx * dx;
        new_cov[1] += dx * dy;
        new_cov[2] += dx * dz;
        new_cov[3] += dy * dy;
        new_cov[4] += dy * dz;
        new_cov[5] += dz * dz;
      }
      for (double& c : new_cov) {
        c /= denom;
      }
    }

    const double shift = std::sqrt(
        (new_center[0] - center[0]) * (new_center[0] - center[0]) +
        (new_center[1] - center[1]) * (new_center[1] - center[1]) +
        (new_center[2] - center[2]) * (new_center[2] - center[2]));
    const bool mask_unchanged = (iter > 0U && next_active == active);

    active = std::move(next_active);
    center = new_center;
    cov = new_cov;
    kept = next_kept;

    scale[0] = std::max(std::sqrt(std::max(0.0, cov[0])),
                        static_cast<double>(cfg.min_scale_cm[0]));
    scale[1] = std::max(std::sqrt(std::max(0.0, cov[3])),
                        static_cast<double>(cfg.min_scale_cm[1]));
    scale[2] = std::max(std::sqrt(std::max(0.0, cov[5])),
                        static_cast<double>(cfg.min_scale_cm[2]));

    if (mask_unchanged && shift <= cfg.convergence_tol_cm) {
      result.converged = true;
      break;
    }
  }

  result.accepted_count = kept;
  for (std::size_t k = 0; k < 3; ++k) {
    result.position_cm[k] = static_cast<float>(center[k]);
  }
  for (std::size_t i = 0; i < 6; ++i) {
    result.covariance_cm2[i] = static_cast<float>(cov[i]);
  }
  if (kept > 0U) {
    const double sqrt_n = std::sqrt(static_cast<double>(kept));
    const double sx = std::sqrt(std::max(0.0, cov[0]));
    const double sy = std::sqrt(std::max(0.0, cov[3]));
    const double sz = std::sqrt(std::max(0.0, cov[5]));
    result.beam_size_cm = {static_cast<float>(sx),
                           static_cast<float>(sy),
                           static_cast<float>(sz)};
    result.uncertainty_cm = {static_cast<float>(sx / sqrt_n),
                             static_cast<float>(sy / sqrt_n),
                             static_cast<float>(sz / sqrt_n)};
  }

  if (kept >= 5U && result.converged) {
    result.quality_flag = 0;
  } else if (kept >= 2U) {
    result.quality_flag = 1;
  } else {
    result.quality_flag = 2;
  }
  return result;
}

std::optional<float> apply_vd_pedestal_gain(float raw_adc,
                                            const VdChannelCalibration& calib) noexcept {
  if (!std::isfinite(raw_adc) || !std::isfinite(calib.pedestal) ||
      !std::isfinite(calib.relative_gain) || calib.relative_gain <= 0.0f) {
    return std::nullopt;
  }
  if (is_masked_channel(calib.status)) {
    return std::nullopt;
  }
  return (raw_adc - calib.pedestal) / calib.relative_gain;
}

std::vector<VdCalibratedStrip> apply_vd_pedestal_gain(
    std::span<const VdRawStripInput> raw_strips,
    std::span<const VdChannelCalibration> calib_table) {
  std::vector<VdCalibratedStrip> out;
  out.reserve(raw_strips.size());

  for (const auto& raw : raw_strips) {
    const VdChannelCalibration* matched = nullptr;
    for (const auto& entry : calib_table) {
      if (entry.layer == raw.layer &&
          entry.module == raw.module &&
          entry.strip == raw.strip) {
        matched = &entry;
        break;
      }
    }

    VdChannelCalibration effective;
    if (matched != nullptr) {
      effective = *matched;
    } else {
      effective.layer = raw.layer;
      effective.module = raw.module;
      effective.strip = raw.strip;
    }

    const auto corrected = apply_vd_pedestal_gain(raw.raw_adc, effective);
    if (!corrected.has_value()) {
      continue;
    }
    const float safe_noise = (std::isfinite(effective.noise) && effective.noise > 0.0f)
                                 ? effective.noise
                                 : 1.0f;
    out.push_back(VdCalibratedStrip{
        .layer = raw.layer,
        .module = raw.module,
        .strip = raw.strip,
        .pulse_height = *corrected,
        .noise = safe_noise,
        .sn_ratio = *corrected / safe_noise,
        .is_z_side = raw.is_z_side,
    });
  }
  return out;
}

}  // namespace delphipp::reco::calibration
