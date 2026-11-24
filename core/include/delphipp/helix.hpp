// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// DELPHI 5-parameter perigee helix representation and analytical propagation
// to cylindrical radii (R) and endcap planes (z) in a solenoidal field.

#pragma once

#include <optional>
#include "delphipp/geometry.hpp"
#include "delphipp/magnetic_field.hpp"

namespace delphipp {

// Standard DELPHI perigee track parameters (matching TKR / EXX conventions):
//   eps    : signed impact parameter in r-phi [cm]
//   z0     : z coordinate at perigee [cm]
//   theta  : polar angle in (0, pi) [rad]
//   phi0   : azimuthal direction at perigee [rad]
//   inv_r  : signed curvature 1/R [cm^-1] (positive for q*Bz > 0)
struct PerigeeTrack {
  float eps = 0.0f;
  float z0 = 0.0f;
  float theta = 1.5707963f;
  float phi0 = 0.0f;
  float inv_r = 0.0f;

  // Transverse momentum p_T [GeV/c] in field bz_gev_per_cm
  [[nodiscard]] float pt_gev(float bz_gev_per_cm) const noexcept;
  // Total momentum |p| [GeV/c]
  [[nodiscard]] float p_gev(float bz_gev_per_cm) const noexcept;
};

struct ExtrapolatedPoint {
  Vec3 xyz{};          // Position [cm]
  Vec3 dir{};          // Unit momentum direction
  float arc_xy_cm = 0.0f;
};

// Evaluate position and unit direction at transverse arc length s_xy [cm].
[[nodiscard]] ExtrapolatedPoint propagate_arc(const PerigeeTrack& trk,
                                              float s_xy_cm) noexcept;

// Intersect helix with cylinder of radius target_r_cm (outgoing branch s_xy >= 0).
[[nodiscard]] std::optional<ExtrapolatedPoint> propagate_to_cylinder(
    const PerigeeTrack& trk, float target_r_cm) noexcept;

// Intersect helix with plane z = target_z_cm (outgoing branch s_xy >= 0).
[[nodiscard]] std::optional<ExtrapolatedPoint> propagate_to_z(
    const PerigeeTrack& trk, float target_z_cm) noexcept;

}  // namespace delphipp
