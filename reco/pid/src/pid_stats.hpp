// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Internal statistical helpers shared across DELPHI++ PID and b-tagging modules
// (chi2 survival probabilities, inverse chi2 quantiles, and Fisher combination).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>

namespace delphipp::reco::pid::detail {

// Upper-tail probability P(chi2 >= x) for 1 degree of freedom: erfc(sqrt(x / 2)).
[[nodiscard]] inline double chi2_prob_1dof(double x) noexcept {
  if (!(x > 0.0)) {
    return 1.0;
  }
  const double z = std::sqrt(0.5 * x);
  return std::clamp(std::erfc(z), 0.0, 1.0);
}

// Regularized upper incomplete gamma Q(a, x) = Gamma(a, x) / Gamma(a) for a > 0, x >= 0.
// Matches CERNLIB PROB(chi2, ndof) with a = ndof / 2.0, x = chi2 / 2.0.
[[nodiscard]] inline double chi2_prob(double chi2, int ndof) noexcept {
  if (ndof <= 0) {
    return 0.0;
  }
  if (!(chi2 > 0.0)) {
    return 1.0;
  }
  if (ndof == 1) {
    return chi2_prob_1dof(chi2);
  }
  if (ndof == 2) {
    return std::clamp(std::exp(-0.5 * chi2), 0.0, 1.0);
  }

  const double a = 0.5 * static_cast<double>(ndof);
  const double x = 0.5 * chi2;
  const double lgamma_a = std::lgamma(a);

  if (x < a + 1.0) {
    // Series expansion for lower incomplete gamma P(a, x), then Q = 1 - P.
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int iter = 0; iter < 200; ++iter) {
      ap += 1.0;
      del *= x / ap;
      sum += del;
      if (std::fabs(del) <= std::fabs(sum) * 1e-14) {
        break;
      }
    }
    const double log_pref = a * std::log(x) - x - lgamma_a;
    const double p_lower = (log_pref > -700.0) ? sum * std::exp(log_pref) : 0.0;
    return std::clamp(1.0 - p_lower, 0.0, 1.0);
  }

  // Modified Lentz continued fraction for Q(a, x).
  constexpr double kTiny = 1e-300;
  double b = x + 1.0 - a;
  double c = 1.0 / kTiny;
  double d = (std::fabs(b) < kTiny) ? (1.0 / kTiny) : (1.0 / b);
  double h = d;
  for (int i = 1; i <= 200; ++i) {
    const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
    b += 2.0;
    d = an * d + b;
    if (std::fabs(d) < kTiny) {
      d = kTiny;
    }
    c = b + an / c;
    if (std::fabs(c) < kTiny) {
      c = kTiny;
    }
    d = 1.0 / d;
    const double delta = d * c;
    h *= delta;
    if (std::fabs(delta - 1.0) <= 1e-14) {
      break;
    }
  }
  const double log_pref = a * std::log(x) - x - lgamma_a;
  const double q_upper = (log_pref > -700.0) ? std::exp(log_pref) * h : 0.0;
  return std::clamp(q_upper, 0.0, 1.0);
}

// Inverse of chi2_prob_1dof(x) = q, equivalent to CERNLIB DCHINV(q, 1).
// Uses Acklam's rational approximation to Phi^-1(1 - q/2) plus two Halley/Newton steps.
[[nodiscard]] inline double inverse_chi2_1dof(double q) noexcept {
  if (q >= 1.0) {
    return 0.0;
  }
  if (q <= 0.0) {
    return 150.0;
  }
  const double p = std::clamp(0.5 * q, 1e-300, 0.5);

  // Acklam inverse normal CDF for lower tail p in (0, 0.5]
  constexpr double a1 = -3.969683028665376e+01;
  constexpr double a2 =  2.209460984245205e+02;
  constexpr double a3 = -2.759285104469687e+02;
  constexpr double a4 =  1.383577518672690e+02;
  constexpr double a5 = -3.066479806614716e+01;
  constexpr double a6 =  2.506628277459239e+00;

  constexpr double b1 = -5.447609879822406e+01;
  constexpr double b2 =  1.615858368580409e+02;
  constexpr double b3 = -1.556989798598866e+02;
  constexpr double b4 =  6.680131188771972e+01;
  constexpr double b5 = -1.328068155288572e+01;

  constexpr double c1 = -7.784894002430293e-03;
  constexpr double c2 = -3.223964580411365e-01;
  constexpr double c3 = -2.400758277161838e+00;
  constexpr double c4 = -2.549732539343734e+00;
  constexpr double c5 =  4.374664141464968e+00;
  constexpr double c6 =  2.938163982698783e+00;

  constexpr double d1 =  7.784695709041462e-03;
  constexpr double d2 =  3.224671290700398e-01;
  constexpr double d3 =  2.445134137142996e+00;
  constexpr double d4 =  3.754408661907416e+00;

  double z = 0.0;
  if (p < 0.02425) {
    const double t = std::sqrt(-2.0 * std::log(p));
    z = -(((((c1 * t + c2) * t + c3) * t + c4) * t + c5) * t + c6) /
         ((((d1 * t + d2) * t + d3) * t + d4) * t + 1.0);
  } else {
    const double t = p - 0.5;
    const double r = t * t;
    z = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * t /
        (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
  }

  // Two Newton refinements on erfc(|z| / sqrt(2)) = q
  double az = std::fabs(z);
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * std::numbers::pi);
  for (int iter = 0; iter < 2; ++iter) {
    const double f = std::erfc(az * std::numbers::sqrt1_2) - q;
    const double pdf = 2.0 * inv_sqrt_2pi * std::exp(-0.5 * az * az);
    if (pdf > 1e-300) {
      az = std::max(0.0, az + f / pdf);
    }
  }
  return az * az;
}

// Fisher log-likelihood combination of independent uniform p-values in (0, 1):
//   P_E = Prod(p_i) * sum_{j=0}^{N-1} (-ln Prod(p_i))^j / j!
// Matches DELPHI AABTAGXX SUBROUTINE bbPROB (lines 3609-3659) in log-space.
[[nodiscard]] inline double fisher_combined_prob(
    std::span<const double> p_values) noexcept {
  double plog = 0.0;
  std::size_t count = 0;
  for (const double p : p_values) {
    if (p > 0.0 && p < 1.0) {
      plog -= std::log(p);
      ++count;
    }
  }
  if (count == 0) {
    return 1.0;
  }
  if (count == 1) {
    return std::clamp(std::exp(-plog), 1e-30, 1.0);
  }

  double xflog = 0.0;
  double pflog = 0.0;
  for (std::size_t j = 1; j < count; ++j) {
    xflog += std::log(plog / static_cast<double>(j));
    const double diff = xflog - pflog;
    if (diff > 50.0) {
      pflog = xflog;
    } else if (diff > -700.0) {
      pflog += std::log1p(std::exp(diff));
    }
  }

  const double log_prob = -plog + pflog;
  if (log_prob <= -69.0775527898) {  // ~ ln(1e-30)
    return 1e-30;
  }
  return std::clamp(std::exp(log_prob), 1e-30, 1.0);
}

}  // namespace delphipp::reco::pid::detail
