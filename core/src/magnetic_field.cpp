// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/magnetic_field.hpp"

#include <cmath>

namespace delphipp {

FieldVector MagneticField::evaluate(const Vec3& xyz) const noexcept {
  constexpr float kB00 = 1e-10f;
  constexpr float kZSat0 = 370.0f;
  constexpr float kZSat1 = 570.0f;
  constexpr float kRSat = 9.0f;

  const float r = std::hypot(xyz[0], xyz[1]) + kB00;
  const float za = std::fabs(xyz[2]);

  float br = 0.0f;
  float bz = 0.0f;
  if (r <= kSolenoidRadiusCm && za <= kSolenoidHalfLenCm) {
    bz = kSolenoidNominalT;
  }

  const float br_gev = scale_ * (br + kB00) * kFieldConvGeVPerCmT;
  const float bz_gev = scale_ * (bz + kB00) * kFieldConvGeVPerCmT;

  FieldVector out;
  out.bx = br_gev * xyz[0] / r;
  out.by = br_gev * xyz[1] / r;
  out.bz = bz_gev;

  // Superconducting quadrupole for the SAT (ufield.car lines 466-470)
  if (za >= kZSat0 && za <= kZSat1 && r <= kRSat && csat_ != 0.0f) {
    out.bx = csat_ * xyz[1];
    out.by = csat_ * xyz[0];
    out.bz = 0.0f;
  }

  out.bmag = std::sqrt(out.bx * out.bx + out.by * out.by + out.bz * out.bz);
  return out;
}

}  // namespace delphipp
