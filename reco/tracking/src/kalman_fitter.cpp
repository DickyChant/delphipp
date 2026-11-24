// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/kalman_fitter.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace delphipp::reco::tracking {

namespace {

double wrap_angle(double d) noexcept {
  constexpr double kPi = std::numbers::pi_v<double>;
  while (d > kPi) d -= 2.0 * kPi;
  while (d < -kPi) d += 2.0 * kPi;
  return d;
}

// Single scalar Kalman update: z = h * x + v, Var(v) = R
// Uses Joseph form P = (I - K h) P (I - K h)^T + K R K^T
double scalar_update(std::array<double, 5>& x,
                     std::array<double, 25>& P,
                     const std::array<double, 5>& h,
                     double residual,
                     double r_var) noexcept {
  std::array<double, 5> Ph{};
  double s = r_var;
  for (std::size_t i = 0; i < 5; ++i) {
    double acc = 0.0;
    for (std::size_t j = 0; j < 5; ++j) {
      acc += P[i * 5 + j] * h[j];
    }
    Ph[i] = acc;
    s += h[i] * acc;
  }
  if (s <= 1e-18) {
    return 0.0;
  }

  std::array<double, 5> K{};
  for (std::size_t i = 0; i < 5; ++i) {
    K[i] = Ph[i] / s;
    x[i] += K[i] * residual;
  }

  // I - K h
  std::array<double, 25> IKh{};
  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 5; ++j) {
      IKh[i * 5 + j] = (i == j ? 1.0 : 0.0) - K[i] * h[j];
    }
  }

  // tmp = (I - K h) * P
  std::array<double, 25> tmp{};
  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 5; ++j) {
      double acc = 0.0;
      for (std::size_t k = 0; k < 5; ++k) {
        acc += IKh[i * 5 + k] * P[k * 5 + j];
      }
      tmp[i * 5 + j] = acc;
    }
  }

  // P_new = tmp * (I - K h)^T + K * R * K^T
  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 5; ++j) {
      double acc = K[i] * r_var * K[j];
      for (std::size_t k = 0; k < 5; ++k) {
        acc += tmp[i * 5 + k] * IKh[j * 5 + k];
      }
      P[i * 5 + j] = acc;
    }
  }

  return (residual * residual) / s;
}

}  // namespace

double highland_theta0(double p_gev, double x_over_x0) noexcept {
  if (p_gev <= 1e-4 || x_over_x0 <= 1e-9) {
    return 0.0;
  }
  const double sq = std::sqrt(x_over_x0);
  const double corr = std::max(0.1, 1.0 + 0.038 * std::log(x_over_x0));
  return (0.0136 / p_gev) * sq * corr;
}

KalmanFitResult kalman_refit(
    const PerigeeTrack& seed,
    float bz_gev_per_cm,
    std::span<const CylinderMeasurement> measurements) noexcept {
  std::array<double, 5> x = {
      seed.eps, seed.z0, seed.theta, seed.phi0, seed.inv_r,
  };

  // Loose initial prior covariance on (eps, z0, theta, phi0, inv_r)
  std::array<double, 25> P{};
  P[0 * 5 + 0] = 1.0;       // 1 cm^2
  P[1 * 5 + 1] = 25.0;      // 25 cm^2
  P[2 * 5 + 2] = 0.04;      // (0.2 rad)^2
  P[3 * 5 + 3] = 0.04;      // (0.2 rad)^2
  P[4 * 5 + 4] = 1e-4;      // (0.01 cm^-1)^2

  KalmanFitResult res;
  std::int32_t n_meas = 0;

  for (const auto& m : measurements) {
    PerigeeTrack cur{
        .eps = static_cast<float>(x[0]),
        .z0 = static_cast<float>(x[1]),
        .theta = static_cast<float>(x[2]),
        .phi0 = static_cast<float>(x[3]),
        .inv_r = static_cast<float>(x[4]),
    };

    const auto pred = propagate_to_cylinder(cur, m.r_cm);
    if (!pred.has_value()) {
      continue;
    }

    // Add multiple scattering process noise to (theta, phi0)
    const double p_gev = std::min(1000.0, static_cast<double>(cur.p_gev(bz_gev_per_cm)));
    const double th0 = highland_theta0(p_gev, static_cast<double>(m.material_x0));
    const double sin_th = std::max(0.1, std::sin(x[2]));
    P[2 * 5 + 2] += th0 * th0;
    P[3 * 5 + 3] += (th0 * th0) / (sin_th * sin_th);

    // 1. R-phi measurement update: r * phi_pred vs r * phi_meas
    // Linearized derivatives at cylinder radius R:
    //   d(R*phi)/deps ~ 1, d(R*phi)/dphi0 ~ R, d(R*phi)/d(inv_r) ~ 0.5 * R^2
    const double pred_phi = std::atan2(static_cast<double>(pred->xyz[1]),
                                       static_cast<double>(pred->xyz[0]));
    const double r_d = static_cast<double>(m.r_cm);
    const double res_rphi = r_d * wrap_angle(static_cast<double>(m.phi_rad) - pred_phi);
    const std::array<double, 5> h_rphi = {1.0, 0.0, 0.0, r_d, 0.5 * r_d * r_d};
    res.chi2 += scalar_update(x, P, h_rphi, res_rphi,
                              static_cast<double>(m.sigma_rphi_cm) * m.sigma_rphi_cm);
    ++n_meas;

    // 2. Z measurement update: z_pred = z0 + s_xy * cot(theta)
    if (m.has_z) {
      const double s_xy = static_cast<double>(pred->arc_xy_cm);
      const double res_z = static_cast<double>(m.z_cm) - static_cast<double>(pred->xyz[2]);
      const double dz_dtheta = -s_xy / (sin_th * sin_th);
      const std::array<double, 5> h_z = {0.0, 1.0, dz_dtheta, 0.0, 0.0};
      res.chi2 += scalar_update(x, P, h_z, res_z,
                                static_cast<double>(m.sigma_z_cm) * m.sigma_z_cm);
      ++n_meas;
    }

    ++res.hits_used;
  }

  res.track = PerigeeTrack{
      .eps = static_cast<float>(x[0]),
      .z0 = static_cast<float>(x[1]),
      .theta = static_cast<float>(x[2]),
      .phi0 = static_cast<float>(x[3]),
      .inv_r = static_cast<float>(x[4]),
  };
  res.covariance = P;
  res.ndf = std::max<std::int32_t>(0, n_meas - 5);
  return res;
}

}  // namespace delphipp::reco::tracking
