// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/pid/hadron_id.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace delphipp::reco::pid {

namespace {

// BBPULL empirical CDF tables from simana_v94d/src/car/dedxid.car (lines 1535-1619)
inline constexpr std::array<float, 100> kPullX = {
    -3.955000f, -3.865000f, -3.775000f, -3.685000f, -3.595000f,
    -3.505000f, -3.415000f, -3.325000f, -3.235000f, -3.145000f,
    -3.055000f, -2.965000f, -2.875000f, -2.785000f, -2.695000f,
    -2.605000f, -2.515000f, -2.425000f, -2.335000f, -2.245000f,
    -2.155000f, -2.065000f, -1.975000f, -1.885000f, -1.795000f,
    -1.705000f, -1.615000f, -1.525000f, -1.435000f, -1.345000f,
    -1.255000f, -1.165000f, -1.075000f, -0.985000f, -0.895000f,
    -0.805000f, -0.715000f, -0.625000f, -0.535000f, -0.445000f,
    -0.355000f, -0.265000f, -0.175000f, -0.085000f,  0.005000f,
     0.095000f,  0.185000f,  0.275000f,  0.365000f,  0.455000f,
     0.545000f,  0.635000f,  0.725000f,  0.815000f,  0.905000f,
     0.995000f,  1.085000f,  1.175000f,  1.265000f,  1.355000f,
     1.445000f,  1.535000f,  1.625000f,  1.715000f,  1.805000f,
     1.895000f,  1.985000f,  2.075000f,  2.165000f,  2.255000f,
     2.345000f,  2.435000f,  2.525000f,  2.615000f,  2.705000f,
     2.795000f,  2.885000f,  2.975000f,  3.065000f,  3.155000f,
     3.245000f,  3.335000f,  3.425000f,  3.515000f,  3.605000f,
     3.695000f,  3.785000f,  3.875000f,  3.965000f,  4.055000f,
     4.145000f,  4.235001f,  4.325000f,  4.415000f,  4.505000f,
     4.595000f,  4.685000f,  4.775001f,  4.865001f,  4.955001f,
};

inline constexpr std::array<float, 100> kPullYiData = {
    0.000020f, 0.000025f, 0.000035f, 0.000084f, 0.000133f,
    0.000212f, 0.000301f, 0.000390f, 0.000498f, 0.000641f,
    0.000868f, 0.001135f, 0.001514f, 0.001904f, 0.002516f,
    0.003280f, 0.004292f, 0.005678f, 0.007296f, 0.009299f,
    0.011780f, 0.014883f, 0.019042f, 0.023812f, 0.029677f,
    0.036426f, 0.045088f, 0.055033f, 0.066616f, 0.080350f,
    0.096037f, 0.113840f, 0.133933f, 0.156033f, 0.180264f,
    0.207228f, 0.236155f, 0.266923f, 0.299836f, 0.333805f,
    0.369387f, 0.406124f, 0.443581f, 0.481240f, 0.518899f,
    0.555635f, 0.592056f, 0.626281f, 0.660803f, 0.693608f,
    0.724656f, 0.753461f, 0.780696f, 0.805810f, 0.828576f,
    0.849650f, 0.868569f, 0.885889f, 0.901003f, 0.914653f,
    0.926714f, 0.937276f, 0.946619f, 0.954473f, 0.961014f,
    0.966751f, 0.971842f, 0.976282f, 0.980006f, 0.983247f,
    0.985980f, 0.988161f, 0.990070f, 0.991500f, 0.992645f,
    0.993691f, 0.994569f, 0.995269f, 0.995901f, 0.996473f,
    0.996907f, 0.997307f, 0.997603f, 0.997894f, 0.998165f,
    0.998407f, 0.998589f, 0.998767f, 0.998925f, 0.999068f,
    0.999196f, 0.999339f, 0.999438f, 0.999571f, 0.999640f,
    0.999699f, 0.999758f, 0.999852f, 0.999936f, 1.000000f,
};

inline constexpr std::array<float, 100> kPullYiMc = {
    0.000013f, 0.000030f, 0.000038f, 0.000063f, 0.000097f,
    0.000139f, 0.000224f, 0.000325f, 0.000456f, 0.000595f,
    0.000764f, 0.001021f, 0.001443f, 0.001975f, 0.002637f,
    0.003625f, 0.004908f, 0.006478f, 0.008545f, 0.011107f,
    0.014457f, 0.018466f, 0.023640f, 0.029518f, 0.036637f,
    0.045394f, 0.055868f, 0.068046f, 0.081909f, 0.097295f,
    0.115090f, 0.134654f, 0.156711f, 0.181343f, 0.207523f,
    0.235759f, 0.265695f, 0.297192f, 0.330045f, 0.363901f,
    0.398918f, 0.434167f, 0.469771f, 0.505416f, 0.541197f,
    0.576788f, 0.610948f, 0.644547f, 0.676526f, 0.707483f,
    0.735862f, 0.762528f, 0.787459f, 0.810909f, 0.832473f,
    0.852341f, 0.870136f, 0.886400f, 0.900608f, 0.913175f,
    0.925088f, 0.935211f, 0.943757f, 0.951471f, 0.958189f,
    0.963852f, 0.968954f, 0.973296f, 0.977094f, 0.980327f,
    0.983179f, 0.985610f, 0.987783f, 0.989589f, 0.991100f,
    0.992337f, 0.993522f, 0.994514f, 0.995189f, 0.995843f,
    0.996422f, 0.996953f, 0.997409f, 0.997747f, 0.998012f,
    0.998232f, 0.998477f, 0.998683f, 0.998865f, 0.999017f,
    0.999173f, 0.999312f, 0.999405f, 0.999532f, 0.999637f,
    0.999730f, 0.999793f, 0.999835f, 0.999937f, 1.000000f,
};

inline constexpr std::array<float, 100> kPullYfGauss = {
    0.000038f, 0.000056f, 0.000080f, 0.000114f, 0.000162f,
    0.000228f, 0.000319f, 0.000442f, 0.000608f, 0.000830f,
    0.001125f, 0.001513f, 0.002020f, 0.002676f, 0.003519f,
    0.004594f, 0.005952f, 0.007654f, 0.009772f, 0.012384f,
    0.015581f, 0.019461f, 0.024134f, 0.029715f, 0.036327f,
    0.044097f, 0.053155f, 0.063629f, 0.075643f, 0.089312f,
    0.104739f, 0.122009f, 0.141187f, 0.162312f, 0.185393f,
    0.210409f, 0.237304f, 0.265984f, 0.296324f, 0.328160f,
    0.361295f, 0.395504f, 0.430539f, 0.466130f, 0.501994f,
    0.537842f, 0.573385f, 0.608341f, 0.642443f, 0.675445f,
    0.707123f, 0.737286f, 0.765774f, 0.792463f, 0.817267f,
    0.840132f, 0.861039f, 0.880002f, 0.897064f, 0.912291f,
    0.925771f, 0.937608f, 0.947919f, 0.956827f, 0.964463f,
    0.970954f, 0.976428f, 0.981007f, 0.984806f, 0.987933f,
    0.990486f, 0.992554f, 0.994215f, 0.995539f, 0.996585f,
    0.997405f, 0.998043f, 0.998535f, 0.998912f, 0.999197f,
    0.999413f, 0.999574f, 0.999693f, 0.999780f, 0.999844f,
    0.999890f, 0.999923f, 0.999947f, 0.999963f, 0.999975f,
    0.999983f, 0.999989f, 0.999992f, 0.999995f, 0.999997f,
    0.999998f, 0.999999f, 0.999999f, 0.999999f, 1.000000f,
};

[[nodiscard]] std::uint8_t compute_hadron_tag(float p_sig_vs_pi,
                                              float p_sig_vs_other,
                                              bool is_best) noexcept {
  if (p_sig_vs_pi > 0.90f && p_sig_vs_other > 0.80f && is_best) {
    return 3;  // Tight
  }
  if (p_sig_vs_pi > 0.70f && p_sig_vs_other > 0.60f) {
    return 2;  // Standard
  }
  if (p_sig_vs_pi > 0.50f) {
    return 1;  // Loose
  }
  return 0;    // None
}

}  // namespace

DedxExpectation bethe_bloch_dedx(float momentum_gev, float mass_gev) noexcept {
  if (!(momentum_gev > 0.001f) || !(mass_gev > 0.0f)) {
    return DedxExpectation{};
  }

  // Exact GETDNEW parameters from dedxid.car (lines 4125-4183)
  constexpr double par[8] = {
      0.13914, 0.65853, 0.93177, 1.00000,
      0.01500, 5.83500, 1.59530, 1.41660,
  };
  constexpr double kNorm = 0.921605;
  constexpr double kN = 3.0;

  const double p = static_cast<double>(momentum_gev);
  const double m = static_cast<double>(mass_gev);
  double x = p / m;
  if (x > 1.0e5) {
    x = 1.0e5;
  }

  const double a1 = par[4];
  const double a2 = par[3];
  const double xa = par[5];
  const double ratio = a2 / (kN * a1);
  const double x0 = xa - (a1 / a2) * std::pow(ratio, kN / (kN - 1.0));
  const double x1 = x0 + std::pow(ratio, 1.0 / (kN - 1.0));

  const double y = std::log(x);
  double add = 0.0;
  double der_add = 0.0;
  if (y < x0) {
    add = 0.0;
    der_add = 0.0;
  } else if (y < x1) {
    const double dx = x1 - y;
    add = a2 * (y - xa) + a1 * std::pow(dx, kN);
    der_add = (a2 - a1 * kN * std::pow(dx, kN - 1.0)) / x;
  } else {
    add = a2 * (y - xa);
    der_add = a2 / x;
  }

  constexpr double am = 2.0;
  const double an = par[6];
  const double bet = par[7];
  const double xx0 = std::pow((am * par[2]) / (an * bet), (am + 1.0) / (an + 1.0));
  const double bet1 = par[2] / (xx0 * xx0) + par[1] - bet / std::pow(xx0, an);

  double rest = 0.0;
  double der_rest = 0.0;
  if (x > xx0) {
    rest = par[1] + par[2] / (x * x);
    der_rest = -par[2] * am * std::pow(x, -am - 1.0);
  } else {
    rest = bet1 + bet / std::pow(x, an);
    der_rest = -bet * an * std::pow(x, -an - 1.0);
  }

  const double dedx =
      ((par[0] * (1.0 + x * x) * (std::log(x) - add)) / (x * x) + rest) / kNorm;
  const double dsurdp_x =
      -(2.0 * par[0] / (x * x * x)) * (std::log(x) - add) +
      par[0] * (1.0 + 1.0 / (x * x)) * (1.0 / x - der_add) +
      der_rest;
  const double dsurdp = dsurdp_x / (m * kNorm);

  return DedxExpectation{
      .dedx_mip = static_cast<float>(dedx),
      .ddedx_dp = static_cast<float>(dsurdp),
  };
}

DedxExpectation bethe_bloch_dedx(float momentum_gev,
                                 rich::Species species) noexcept {
  const auto idx = static_cast<std::size_t>(species);
  if (idx >= rich::kSpeciesMassGeV.size()) {
    return DedxExpectation{};
  }
  return bethe_bloch_dedx(momentum_gev, rich::kSpeciesMassGeV[idx]);
}

float expected_dedx_error(float expected_dedx_mip,
                          std::int32_t num_wires,
                          float gap_cm,
                          bool edge_region,
                          std::int32_t* quality_out) noexcept {
  if (quality_out != nullptr) {
    *quality_out = -1;
  }
  if (num_wires < 1 || !(expected_dedx_mip > 0.0f)) {
    return 1.0f;
  }

  const double nw = static_cast<double>(num_wires);
  const double dedx_th = static_cast<double>(expected_dedx_mip);
  const double gr = std::clamp(static_cast<double>(gap_cm) / 0.4, 0.1, 1.71);

  double error =
      dedx_th *
      std::sqrt((1.0757 / nw + 0.823 / (nw * nw)) /
                std::pow(dedx_th, 1.29) / std::pow(gr, 0.926)) *
      1.0286;

  if (edge_region) {
    error *= 1.33;
  }

  std::int32_t quality = 0;
  if (num_wires > 37) {
    const double rel = error / dedx_th;
    if (rel < 0.06) {
      quality = 4;
    } else if (rel < 0.07) {
      quality = 3;
    } else if (rel < 0.09) {
      quality = 2;
    } else if (rel < 0.11) {
      quality = 1;
    }
  }
  if (quality_out != nullptr) {
    *quality_out = quality;
  }

  return static_cast<float>(error);
}

float gaussianize_dedx_pull(float raw_pull, bool is_simulation) noexcept {
  if (!std::isfinite(raw_pull)) {
    return raw_pull;
  }
  if (raw_pull <= kPullX.front() || raw_pull >= kPullX.back()) {
    return raw_pull;
  }

  const auto& yi = is_simulation ? kPullYiMc : kPullYiData;

  std::size_t i1 = 0;
  std::size_t i2 = 0;
  for (std::size_t i = 0; i < kPullX.size(); ++i) {
    if (raw_pull > kPullX[i]) {
      i1 = i;
    }
    if (raw_pull <= kPullX[i]) {
      i2 = i;
      break;
    }
  }

  float y = yi[i1];
  if (i1 != i2) {
    const float dx = kPullX[i1] - kPullX[i2];
    if (std::fabs(dx) > 1e-12f) {
      y = ((yi[i1] - yi[i2]) * raw_pull + yi[i2] * kPullX[i1] -
           yi[i1] * kPullX[i2]) /
          dx;
    }
  }

  std::size_t j1 = 0;
  std::size_t j2 = kPullYfGauss.size() - 1;
  bool found = false;
  for (std::size_t i = 0; i < kPullYfGauss.size(); ++i) {
    if (y > kPullYfGauss[i]) {
      j1 = i;
    }
    if (y <= kPullYfGauss[i]) {
      j2 = i;
      found = true;
      break;
    }
  }
  if (!found) {
    j2 = kPullYfGauss.size() - 1;
  }

  if (j1 == j2) {
    return kPullX[j1];
  }
  const float dy = kPullYfGauss[j1] - kPullYfGauss[j2];
  if (dy >= 0.0f) {
    return raw_pull;
  }
  return (y * (kPullX[j1] - kPullX[j2]) -
          (kPullYfGauss[j2] * kPullX[j1] - kPullYfGauss[j1] * kPullX[j2])) /
         dy;
}

float dedx_pairwise_probability(float momentum_gev,
                                float measured_dedx,
                                std::int32_t num_wires,
                                rich::Species signal_species,
                                rich::Species background_species,
                                float gap_cm,
                                bool edge_region) noexcept {
  if (!(measured_dedx > 0.0f) || num_wires < 1) {
    return -0.5f;
  }
  const auto exp1 = bethe_bloch_dedx(momentum_gev, signal_species);
  const auto exp2 = bethe_bloch_dedx(momentum_gev, background_species);
  if (!(exp1.dedx_mip > 0.0f) || !(exp2.dedx_mip > 0.0f)) {
    return -0.5f;
  }

  const float err1 =
      expected_dedx_error(exp1.dedx_mip, num_wires, gap_cm, edge_region);
  const float err2 =
      expected_dedx_error(exp2.dedx_mip, num_wires, gap_cm, edge_region);
  const float pull1 =
      gaussianize_dedx_pull((measured_dedx - exp1.dedx_mip) / err1);
  const float pull2 =
      gaussianize_dedx_pull((measured_dedx - exp2.dedx_mip) / err2);

  if (std::fabs(pull1) > 5.0f) {
    return 0.0f;
  }
  if (std::fabs(pull2) > 5.0f && (std::fabs(pull2) - std::fabs(pull1)) > 2.0f) {
    return 1.0f;
  }
  const double pb1 = std::exp(-0.5 * static_cast<double>(pull1) * pull1);
  const double pb2 = std::exp(-0.5 * static_cast<double>(pull2) * pull2);
  const double denom = pb1 + pb2;
  if (!(denom > 0.0)) {
    return 0.0f;
  }
  return static_cast<float>(pb1 / denom);
}

HadronIdResult identify_hadron(const HadronIdInput& input) noexcept {
  HadronIdResult res{};

  const bool has_dedx = (input.tpc_dedx > 0.0f && input.tpc_num_wires > 0);
  std::optional<rich::RichLikelihoodResult> liq_res;
  std::optional<rich::RichLikelihoodResult> gas_res;

  if (input.rich_liquid.has_value()) {
    liq_res = rich::evaluate_hypotheses(*input.rich_liquid);
  }
  if (input.rich_gas.has_value()) {
    gas_res = rich::evaluate_hypotheses(*input.rich_gas);
  }

  double min_chi2 = std::numeric_limits<double>::infinity();
  std::size_t best_idx = static_cast<std::size_t>(rich::Species::Pion);

  for (std::size_t i = 0; i < 5; ++i) {
    const auto sp = static_cast<rich::Species>(i);
    const auto bb = bethe_bloch_dedx(input.momentum_gev, sp);
    res.expected_dedx[i] = bb.dedx_mip;

    if (has_dedx && bb.dedx_mip > 0.0f) {
      std::int32_t q = -1;
      const float err = expected_dedx_error(
          bb.dedx_mip, input.tpc_num_wires, input.tpc_gap_cm,
          input.tpc_edge_region, &q);
      res.dedx_error[i] = err;
      if (i == static_cast<std::size_t>(rich::Species::Pion)) {
        res.dedx_quality = q;
      }
      const float raw_pull = (input.tpc_dedx - bb.dedx_mip) / err;
      const float pull =
          input.gaussianize_dedx_pull
              ? gaussianize_dedx_pull(raw_pull, input.is_simulation)
              : raw_pull;
      res.dedx_pull[i] = pull;
      res.dedx_chi2[i] = pull * pull;
    }

    double r_chi2 = 0.0;
    if (input.rich_chi2_override.has_value()) {
      r_chi2 += static_cast<double>((*input.rich_chi2_override)[i]);
    } else {
      if (liq_res.has_value()) {
        r_chi2 += static_cast<double>(liq_res->chi2[i]);
      }
      if (gas_res.has_value()) {
        r_chi2 += static_cast<double>(gas_res->chi2[i]);
      }
    }
    res.rich_chi2[i] = static_cast<float>(r_chi2);

    const double tot = static_cast<double>(res.dedx_chi2[i]) + r_chi2;
    res.total_chi2[i] = static_cast<float>(tot);
    if (tot < min_chi2) {
      min_chi2 = tot;
      best_idx = i;
    }
  }

  std::array<double, 5> raw_l{};
  double sum_l = 0.0;
  for (std::size_t i = 0; i < 5; ++i) {
    const double d = static_cast<double>(res.total_chi2[i]) - min_chi2;
    raw_l[i] = (d < 500.0) ? std::exp(-0.5 * d) : 0.0;
    sum_l += raw_l[i];
  }
  if (sum_l > 0.0) {
    for (std::size_t i = 0; i < 5; ++i) {
      res.likelihood[i] = static_cast<float>(raw_l[i] / sum_l);
    }
  } else {
    res.likelihood[static_cast<std::size_t>(rich::Species::Pion)] = 1.0f;
  }

  res.best_species = static_cast<rich::Species>(best_idx);

  const float l_pi = res.likelihood[static_cast<std::size_t>(rich::Species::Pion)];
  const float l_k = res.likelihood[static_cast<std::size_t>(rich::Species::Kaon)];
  const float l_p = res.likelihood[static_cast<std::size_t>(rich::Species::Proton)];

  res.prob_kaon_vs_pion = (l_k + l_pi > 0.0f) ? (l_k / (l_k + l_pi)) : 0.5f;
  res.prob_proton_vs_pion = (l_p + l_pi > 0.0f) ? (l_p / (l_p + l_pi)) : 0.5f;
  res.prob_proton_vs_kaon = (l_p + l_k > 0.0f) ? (l_p / (l_p + l_k)) : 0.5f;

  res.kaon_tag = compute_hadron_tag(
      res.prob_kaon_vs_pion, 1.0f - res.prob_proton_vs_kaon,
      res.best_species == rich::Species::Kaon);
  res.proton_tag = compute_hadron_tag(
      res.prob_proton_vs_pion, res.prob_proton_vs_kaon,
      res.best_species == rich::Species::Proton);

  return res;
}

}  // namespace delphipp::reco::pid
