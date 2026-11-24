// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/magnetic_field.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace delphipp {

namespace {

// Linear interpolation matching FINT(X, Y1, X1, Y2, X2) in ufield.car line 334
bool interpolate_1d(float x,
                    const std::vector<float>& grid,
                    const std::vector<float>& values,
                    float& out) noexcept {
  if (grid.size() < 2 || grid.size() != values.size()) return false;
  if (x < grid.front() || x > grid.back()) return false;
  const auto it = std::upper_bound(grid.begin(), grid.end(), x);
  std::size_t i2 = static_cast<std::size_t>(it - grid.begin());
  if (i2 == 0) i2 = 1;
  if (i2 >= grid.size()) i2 = grid.size() - 1;
  const std::size_t i1 = i2 - 1;
  const float dx = grid[i2] - grid[i1];
  if (std::fabs(dx) < 1e-12f) {
    out = values[i1];
  } else {
    out = values[i1] + (values[i2] - values[i1]) * (x - grid[i1]) / dx;
  }
  return true;
}

}  // namespace

void MagneticField::set_barrel_layers(std::vector<FieldLayerTable> layers) {
  barrel_layers_ = std::move(layers);
}

void MagneticField::set_endcap_layers(std::vector<FieldLayerTable> layers) {
  endcap_layers_ = std::move(layers);
}

FieldVector MagneticField::evaluate(const Vec3& xyz) const noexcept {
  constexpr float kB00 = 1e-10f;
  constexpr float kPi = std::numbers::pi_v<float>;
  constexpr float kAlfa1B = kPi / 24.0f;
  constexpr float kAlfa2B = 2.0f * kAlfa1B;
  constexpr float kAlfa1E = kPi / 12.0f;
  constexpr float kAlfa2E = 2.0f * kAlfa1E;
  constexpr float kZSat0 = 370.0f;
  constexpr float kZSat1 = 570.0f;
  constexpr float kRSat = 9.0f;
  constexpr float kTzanko = 1.0f;

  const float r = std::hypot(xyz[0], xyz[1]) + kB00;
  const float za = std::fabs(xyz[2]);

  float br = 0.0f;
  float bz = 0.0f;
  bool found = false;

  // 1. Barrel grid lookup (ufield.car lines 375-413)
  if (!barrel_layers_.empty()) {
    float x1 = xyz[0];
    float y1 = xyz[1];
    if (accspc_) {
      x1 = std::copysign(std::max(0.0001f, std::fabs(x1) - kTzanko), x1);
      y1 = std::copysign(std::max(0.0001f, std::fabs(y1) - kTzanko), y1);
    }
    const float rr = std::hypot(x1, y1) + kB00;
    if (x1 == 0.0f && y1 == 0.0f) x1 = 0.0001f;
    const float phi = std::fabs(std::atan2(y1, x1));
    const float dev = std::fmod(phi, kAlfa2B) - kAlfa1B;
    const float dev2 = dev * dev;
    const float rpseu = rr * (1.0f - dev2 * (0.5f - dev2 / 24.0f));

    for (const auto& lay : barrel_layers_) {
      if (rpseu >= lay.coord_min_cm && rpseu <= lay.coord_max_cm) {
        float bz_t = 0.0f, br_t = 0.0f;
        if (interpolate_1d(za, lay.bz_grid_cm, lay.bz_tesla, bz_t) &&
            interpolate_1d(za, lay.br_grid_cm, lay.br_tesla, br_t)) {
          bz = bz_t;
          br = br_t;
          found = true;
          break;
        }
      }
    }
  }

  // 2. Endcap grid lookup (ufield.car lines 416-445)
  if (!found && !endcap_layers_.empty()) {
    float x1 = xyz[0];
    float y1 = xyz[1];
    if (x1 == 0.0f && y1 == 0.0f) x1 = 0.0001f;
    const float phi = std::fabs(std::atan2(y1, x1)) + kAlfa1E;
    const float dev = std::fmod(phi, kAlfa2E) - kAlfa1E;
    const float dev2 = dev * dev;
    const float rpseu = r * (1.0f - dev2 * (0.5f - dev2 / 24.0f));

    for (const auto& lay : endcap_layers_) {
      if (za >= lay.coord_min_cm && za <= lay.coord_max_cm) {
        float bz_t = 0.0f, br_t = 0.0f;
        if (interpolate_1d(rpseu, lay.bz_grid_cm, lay.bz_tesla, bz_t) &&
            interpolate_1d(rpseu, lay.br_grid_cm, lay.br_tesla, br_t)) {
          bz = bz_t;
          br = br_t;
          found = true;
          break;
        }
      }
    }
  }

  if (found) {
    if (xyz[2] < 0.0f) br = -br;
  } else if (r <= kSolenoidRadiusCm && za <= kSolenoidHalfLenCm) {
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
