// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI magnetic field model (simana_v94d/src/car/ufield.car).
// Returns field in DELPHI units of GeV/cm (B_GeV_per_cm * R_cm = p_T_GeV),
// with FCT = 0.00299793 GeV/(cm*T) and nominal solenoid B0 = 1.2 T.

#pragma once

#include <array>
#include "delphipp/geometry.hpp"

namespace delphipp {

inline constexpr float kFieldConvGeVPerCmT = 0.00299793f;  // FCT in ufield.car
inline constexpr float kSolenoidNominalT = 1.2f;           // B0
inline constexpr float kSolenoidRadiusCm = 275.0f;         // RSOL
inline constexpr float kSolenoidHalfLenCm = 370.0f;        // ZMAXEN nominal

struct FieldVector {
  float bx = 0.0f;  // GeV/cm
  float by = 0.0f;  // GeV/cm
  float bz = 0.0f;  // GeV/cm
  float bmag = 0.0f;

  [[nodiscard]] float bz_tesla() const noexcept { return bz / kFieldConvGeVPerCmT; }
};

class MagneticField {
 public:
  MagneticField() = default;

  // Mirrors UFCSCL: scale nominal current/field (default 1.0)
  void set_scale(float scale) noexcept { scale_ = scale; }
  [[nodiscard]] float scale() const noexcept { return scale_; }

  // Mirrors SAT superconducting quadrupole strength CSAT (GeV/cm^2)
  void set_sat_quad_gradient(float csat) noexcept { csat_ = csat; }

  // Mirrors SUBROUTINE UFIELD(XYZ, BXYZ): evaluates field at DELPHI (x,y,z) [cm]
  [[nodiscard]] FieldVector evaluate(const Vec3& xyz_cm) const noexcept;

 private:
  float scale_ = 1.0f;
  float csat_ = 0.0f;  // SAT quadrupole gradient (370 <= |z| <= 570, r <= 9 cm)
};

}  // namespace delphipp
