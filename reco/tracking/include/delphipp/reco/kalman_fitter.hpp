// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 5-parameter Kalman track refitter with Highland multiple
// scattering process noise (replaces core of tkfb.car / tkff.car / tkrefitp.car).

#pragma once

#include <array>
#include <cstdint>
#include <span>

#include "delphipp/helix.hpp"

namespace delphipp::reco::tracking {

// One cylindrical tracking measurement (from VD, ID, TPC, or OD)
struct CylinderMeasurement {
  float r_cm = 0.0f;
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
  float sigma_rphi_cm = 0.02f;
  float sigma_z_cm = 0.10f;
  bool has_z = true;
  float material_x0 = 0.005f;  // Radiation lengths x/X0 crossed at this layer
};

struct KalmanFitResult {
  PerigeeTrack track{};
  std::array<double, 25> covariance{};  // 5x5 symmetric covariance matrix
  double chi2 = 0.0;
  std::int32_t ndf = 0;
  std::int32_t hits_used = 0;
};

// Highland multiple-scattering RMS angle theta_0 [rad] for momentum p [GeV/c]
[[nodiscard]] double highland_theta0(double p_gev, double x_over_x0) noexcept;

// Refit a seed PerigeeTrack against a sequence of cylindrical measurements
// using sequential Kalman updates in Joseph stabilized form.
[[nodiscard]] KalmanFitResult kalman_refit(
    const PerigeeTrack& seed,
    float bz_gev_per_cm,
    std::span<const CylinderMeasurement> measurements) noexcept;

}  // namespace delphipp::reco::tracking
