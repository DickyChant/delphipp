// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>

#include "delphipp/helix.hpp"
#include "delphipp/magnetic_field.hpp"

int main() {
  using namespace delphipp;

  // 1. MagneticField inside solenoid (B0 = 1.2 T)
  MagneticField field;
  const auto center = field.evaluate({0.0f, 0.0f, 0.0f});
  assert(std::fabs(center.bz_tesla() - 1.2f) < 1e-4f);
  assert(std::fabs(center.bz - 1.2f * kFieldConvGeVPerCmT) < 1e-6f);

  // Outside solenoid
  const auto outside = field.evaluate({300.0f, 0.0f, 0.0f});
  assert(std::fabs(outside.bz_tesla()) < 1e-6f);

  // SAT superconducting quadrupole region
  field.set_sat_quad_gradient(0.01f);
  const auto sat = field.evaluate({3.0f, 4.0f, 450.0f});
  assert(std::fabs(sat.bx - 0.04f) < 1e-6f);
  assert(std::fabs(sat.by - 0.03f) < 1e-6f);
  assert(sat.bz == 0.0f);

  // 2. PerigeeTrack momentum & cylinder propagation
  const float p_t = 10.0f;  // 10 GeV/c transverse momentum
  PerigeeTrack trk{
      .eps = 0.5f,
      .z0 = 2.0f,
      .theta = 1.2f,
      .phi0 = 0.4f,
      .inv_r = center.bz / p_t,
  };
  assert(std::fabs(trk.pt_gev(center.bz) - p_t) < 1e-4f);

  const auto at_tof = propagate_to_cylinder(trk, 311.25f);
  assert(at_tof.has_value());
  const float r_hit = std::hypot(at_tof->xyz[0], at_tof->xyz[1]);
  assert(std::fabs(r_hit - 311.25f) < 1e-2f);
  const float dir_norm = std::sqrt(at_tof->dir[0] * at_tof->dir[0] +
                                   at_tof->dir[1] * at_tof->dir[1] +
                                   at_tof->dir[2] * at_tof->dir[2]);
  assert(std::fabs(dir_norm - 1.0f) < 1e-5f);

  // 3. Endcap plane propagation (z = 450 cm)
  const auto at_hof = propagate_to_z(trk, 450.0f);
  assert(at_hof.has_value());
  assert(std::fabs(at_hof->xyz[2] - 450.0f) < 1e-3f);

  std::cout << "test_field_helix: PASS\n";
  return 0;
}
