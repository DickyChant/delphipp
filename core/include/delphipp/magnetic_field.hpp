// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI magnetic field model (simana_v94d/src/car/ufield.car).
// Supports:
//   - Nominal solenoid (B0 = 1.2 T, RSOL = 275 cm, FCT = 0.00299793 GeV/(cm*T))
//   - Barrel 24-fold & Endcap 12-fold pseudo-radius grid interpolation (IDLAY + FINT)
//   - 2 cm HAC iron-spacer correction (ACCSPC / UFSPAC)
//   - SAT superconducting quadrupole (370 <= |z| <= 570 cm, r <= 9 cm)

#pragma once

#include <array>
#include <vector>
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

// One radial barrel layer or axial endcap layer of the UFIELD lookup table
struct FieldLayerTable {
  float coord_min_cm = 0.0f;
  float coord_max_cm = 0.0f;
  std::vector<float> grid_cm;   // Monotonic |z| (barrel) or r_pseudo (endcap)
  std::vector<float> bz_tesla;  // B_z [T] at grid points
  std::vector<float> br_tesla;  // B_r [T] at grid points
};

class MagneticField {
 public:
  MagneticField() = default;

  void set_scale(float scale) noexcept { scale_ = scale; }
  [[nodiscard]] float scale() const noexcept { return scale_; }

  void set_iron_spacers(bool enabled) noexcept { accspc_ = enabled; }
  void set_sat_quad_gradient(float csat) noexcept { csat_ = csat; }

  // Load optional barrel (indexed by r_pseudo) and endcap (indexed by |z|) field maps
  void set_barrel_layers(std::vector<FieldLayerTable> layers);
  void set_endcap_layers(std::vector<FieldLayerTable> layers);

  // Mirrors SUBROUTINE UFIELD(XYZ, BXYZ): evaluates field at DELPHI (x,y,z) [cm]
  [[nodiscard]] FieldVector evaluate(const Vec3& xyz_cm) const noexcept;

 private:
  float scale_ = 1.0f;
  float csat_ = 0.0f;
  bool accspc_ = false;
  std::vector<FieldLayerTable> barrel_layers_;
  std::vector<FieldLayerTable> endcap_layers_;
};

}  // namespace delphipp
