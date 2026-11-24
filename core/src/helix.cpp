// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/helix.hpp"

#include <cmath>
#include <limits>

namespace delphipp {

float PerigeeTrack::pt_gev(float bz_gev_per_cm) const noexcept {
  const float k = std::fabs(inv_r);
  if (k < 1e-12f) {
    return std::numeric_limits<float>::infinity();
  }
  return std::fabs(bz_gev_per_cm) / k;
}

float PerigeeTrack::p_gev(float bz_gev_per_cm) const noexcept {
  const float st = std::sin(theta);
  if (std::fabs(st) < 1e-9f) {
    return std::numeric_limits<float>::infinity();
  }
  return pt_gev(bz_gev_per_cm) / std::fabs(st);
}

ExtrapolatedPoint propagate_arc(const PerigeeTrack& trk, float s_xy) noexcept {
  const double c0 = std::cos(static_cast<double>(trk.phi0));
  const double s0 = std::sin(static_cast<double>(trk.phi0));
  const double cot_theta = 1.0 / std::tan(static_cast<double>(trk.theta));
  const double sin_theta = std::sin(static_cast<double>(trk.theta));
  const double cos_theta = std::cos(static_cast<double>(trk.theta));

  // Perigee point: (-eps * sin(phi0), eps * cos(phi0), z0)
  const double x0 = -static_cast<double>(trk.eps) * s0;
  const double y0 = static_cast<double>(trk.eps) * c0;
  const double k = static_cast<double>(trk.inv_r);
  const double s = static_cast<double>(s_xy);

  double x = 0.0;
  double y = 0.0;
  double phi = static_cast<double>(trk.phi0) + k * s;

  if (std::fabs(k * s) < 1e-6) {
    // Straight-line limit
    x = x0 + s * c0;
    y = y0 + s * s0;
  } else {
    const double r = 1.0 / k;
    x = x0 + r * (std::sin(phi) - s0);
    y = y0 - r * (std::cos(phi) - c0);
  }

  const double z = static_cast<double>(trk.z0) + s * cot_theta;

  return ExtrapolatedPoint{
      .xyz = {static_cast<float>(x), static_cast<float>(y), static_cast<float>(z)},
      .dir = {static_cast<float>(sin_theta * std::cos(phi)),
              static_cast<float>(sin_theta * std::sin(phi)),
              static_cast<float>(cos_theta)},
      .arc_xy_cm = s_xy,
  };
}

std::optional<ExtrapolatedPoint> propagate_to_cylinder(
    const PerigeeTrack& trk, float target_r_cm) noexcept {
  const double R = static_cast<double>(target_r_cm);
  const double eps = static_cast<double>(trk.eps);
  if (R < std::fabs(eps)) {
    return std::nullopt;
  }

  const double k = static_cast<double>(trk.inv_r);
  double s_xy = 0.0;
  if (std::fabs(k) < 1e-9) {
    s_xy = std::sqrt(R * R - eps * eps);
  } else {
    // Exact chord relation for a helix from perigee (r0 = |eps|, dr/ds = 0):
    //   r^2(s) = eps^2 + (1 - k*eps) * (2/k)^2 * sin^2(k*s/2)
    const double denom = 1.0 - k * eps;
    if (denom <= 0.0) {
      return std::nullopt;
    }
    const double sin2_half = (R * R - eps * eps) * (k * k) / (4.0 * denom);
    if (sin2_half < 0.0 || sin2_half > 1.0) {
      return std::nullopt;
    }
    const double half_angle = std::asin(std::sqrt(sin2_half));
    s_xy = 2.0 * half_angle / std::fabs(k);
  }

  return propagate_arc(trk, static_cast<float>(s_xy));
}

std::optional<ExtrapolatedPoint> propagate_to_z(
    const PerigeeTrack& trk, float target_z_cm) noexcept {
  const double cot_theta = 1.0 / std::tan(static_cast<double>(trk.theta));
  if (std::fabs(cot_theta) < 1e-9) {
    return std::nullopt;
  }
  const double s_xy = (static_cast<double>(target_z_cm) - static_cast<double>(trk.z0)) / cot_theta;
  if (s_xy < 0.0) {
    return std::nullopt;
  }
  return propagate_arc(trk, static_cast<float>(s_xy));
}

}  // namespace delphipp
