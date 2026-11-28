// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/pid/electron_id.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>

#include "delphipp/reco/pid/hadron_id.hpp"
#include "delphipp/reco/rich.hpp"
#include "pid_stats.hpp"

namespace delphipp::reco::pid {

namespace {

// Map 9 HPC readout layers -> 10-layer ELECID indices (layer 4 is trigger gap):
// layers 0..2 -> HPC 1..3 (3 gaps each), layers 3..5 -> HPC 5..7 (4 gaps each),
// layers 6..8 -> HPC 8..10 (6 gaps each).
inline constexpr std::array<int, 9> kHpcGapsPerLayer = {3, 3, 3, 4, 4, 4, 6, 6, 6};

struct ShowerModelParams {
  double v1 = 0.16;
  double v2 = 0.25;
  double energy = 10.0;
};

void compute_expected_9layers(double theta_rad,
                              double z_entry_cm,
                              const ShowerModelParams& p,
                              std::array<double, 9>& out_layers) noexcept {
  const double st = std::max(std::fabs(std::sin(theta_rad)), 0.25);
  const double ct = std::cos(theta_rad) / st;
  const double t0 = 0.7 / st;
  const double dt = 0.430 / st;
  const double dt_cm = 1.20;

  const double v1 = std::max(p.v1, 1e-3);
  const double v2 = std::max(p.v2, 1e-3);
  const double alpha = v2 / std::pow(v1, 1.5) + 1.0;
  const double beta = v2 / std::sqrt(v1);

  const double lgamma_a = std::lgamma(alpha);
  const double pref = std::log(std::max(dt * p.energy * beta, 1e-30)) - lgamma_a;

  double t = t0;
  int sample_idx = 0;
  for (std::size_t layer = 0; layer < 9; ++layer) {
    double acc = 0.0;
    // Account for skipped trigger gap (1 gap between layer 2 and layer 3)
    if (layer == 3) {
      ++sample_idx;
      t += dt;
    }
    const int n_gaps = kHpcGapsPerLayer[layer];
    for (int g = 0; g < n_gaps; ++g) {
      ++sample_idx;
      t += dt;
      // Mild edge-z material attenuation near z=92.4, 183, 247 cm (ELZMAT)
      const double z_abs = std::fabs(z_entry_cm + static_cast<double>(sample_idx) * dt_cm * ct);
      double sensi = 1.0;
      if (std::fabs(z_abs - 92.4) < 1.2 || std::fabs(z_abs - 183.0) < 1.8 ||
          z_abs > 245.0) {
        sensi = 0.7;
      }
      const double flush = (layer == 0) ? 1.2 : 1.0;
      const double xhelp = std::clamp((alpha - 1.0) * std::log(std::max(beta * t, 1e-12)) -
                                          beta * t + pref,
                                      -40.0, 40.0);
      acc += flush * std::exp(xhelp) * sensi;
    }
    out_layers[layer] = acc;
  }
}

struct ShowerFitResult {
  double chi2 = -1.0;
  int ndof = 0;
};

ShowerFitResult fit_hpc_longitudinal_profile(const std::array<float, 9>& measured_layers,
                                             float ecal_energy_gev,
                                             float theta_rad,
                                             float z_entry_cm,
                                             bool is_crack) noexcept {
  // 94D error scale CFACT = 0.90 (elepha32.car line 13862)
  double cfact = 0.90;
  const double z_abs = std::fabs(static_cast<double>(z_entry_cm));
  if (is_crack && z_abs < 200.0) {
    cfact *= 0.9;
  }
  if (z_abs > 83.0 && z_abs < 87.0) {
    cfact *= 0.5;
  }

  const double ee = std::max(static_cast<double>(ecal_energy_gev), 1.0);
  const double log_ee = std::log(ee);

  // Layer energy uncertainties (ELSIGMA, lines 14219-14233)
  constexpr std::array<double, 9> kErrA = {
      0.1150, 0.1752, 0.1386, 0.1294, 0.1611, 0.1766, 0.0596, 0.0371, 0.0161};
  constexpr std::array<double, 9> kErrB = {
      0.3642, 0.4803, 0.6692, 0.7318, 0.5945, 0.5813, 0.9929, 0.9656, 1.0530};
  constexpr std::array<double, 9> kErrDiv = {
      2.1, 2.1, 1.7, 1.7, 1.7, 2.1, 2.1, 3.0, 4.0};

  std::array<double, 9> sigma{};
  std::array<double, 9> elm{};
  for (std::size_t i = 0; i < 9; ++i) {
    const double raw_err = (kErrA[i] * std::pow(ee, kErrB[i])) / kErrDiv[i] / cfact;
    sigma[i] = std::max(raw_err, 0.02);
    elm[i] = std::max(static_cast<double>(measured_layers[i]), 0.0075);
  }

  const double v1_exp = 1.0 / (4.234 + 0.7879 * log_ee);
  const double v2_exp = (1.0 + 0.203e-3 * z_abs) / (2.8453 + 0.3042 * log_ee);
  const double sig_v1 = std::max(1.0 / std::max(-1.408 + 16.730 * log_ee, 1.0) / 1.25 / cfact, 0.01);
  const double sig_v2 = std::max(1.0 / (1.638 + 6.7820 * log_ee) / 1.65 / cfact, 0.01);

  auto eval_chi2 = [&](const ShowerModelParams& p) noexcept -> double {
    std::array<double, 9> pred{};
    compute_expected_9layers(static_cast<double>(theta_rad),
                             static_cast<double>(z_entry_cm), p, pred);
    double sum_chi2 = 0.0;
    double max_chi2 = 0.0;
    for (std::size_t i = 0; i < 9; ++i) {
      const double pull = (elm[i] - pred[i]) / sigma[i];
      const double c = pull * pull;
      sum_chi2 += c;
      if (c > max_chi2) {
        max_chi2 = c;
      }
    }
    // ELFCHIQ drops the single largest-pull layer to suppress dead-channel spikes
    sum_chi2 -= max_chi2;
    const double pv1 = (v1_exp - p.v1) / sig_v1;
    const double pv2 = (v2_exp - p.v2) / sig_v2;
    sum_chi2 += pv1 * pv1 + pv2 * pv2;
    sum_chi2 += 100.0 * std::pow(std::max(p.energy / ee, 2.5) - 2.5, 2.0);
    return sum_chi2;
  };

  // Start at expected electron shape (v1_exp, v2_exp, ee) and refine locally
  ShowerModelParams best{.v1 = v1_exp, .v2 = v2_exp, .energy = ee};
  double best_chi2 = eval_chi2(best);

  // Moment estimate from ELFSTRT as an alternate seed
  const double st = std::max(std::fabs(std::sin(static_cast<double>(theta_rad))), 0.25);
  const double t0 = 0.7 / st;
  const double dt = 0.430 / st;
  constexpr std::array<double, 9> kWir = {2.0, 5.0, 8.0, 12.5, 16.5, 20.5, 25.5, 31.5, 37.5};
  double e0 = 0.0;
  double e1 = 0.0;
  double e2 = 0.0;
  for (std::size_t i = 0; i < 9; ++i) {
    const double t = t0 + kWir[i] * dt;
    e0 += elm[i];
    e1 += elm[i] * t;
    e2 += elm[i] * t * t;
  }
  if (e0 > 0.0) {
    const double tm = e1 / e0;
    const double var = (e2 / e0) - tm * tm;
    if (var > 1e-6) {
      const double b_mom = tm / var;
      const double a_mom = tm * b_mom;
      if (a_mom > 1.05) {
        const ShowerModelParams mom{
            .v1 = b_mom / (a_mom - 1.0),
            .v2 = std::sqrt(std::pow(b_mom, 3.0) / (a_mom - 1.0)),
            .energy = e0,
        };
        const double mom_chi2 = eval_chi2(mom);
        if (mom_chi2 < best_chi2) {
          best = mom;
          best_chi2 = mom_chi2;
        }
      }
    }
  }

  // Deterministic coordinate-descent polish around best seed
  double step_v1 = 0.1 * best.v1;
  double step_v2 = 0.1 * best.v2;
  double step_e = 0.05 * best.energy;
  for (int iter = 0; iter < 16; ++iter) {
    bool improved = false;
    for (int sgn : {-1, 1}) {
      ShowerModelParams cand = best;
      cand.v1 = std::clamp(best.v1 + static_cast<double>(sgn) * step_v1, 1e-3, 1.0);
      const double c = eval_chi2(cand);
      if (c < best_chi2) {
        best = cand;
        best_chi2 = c;
        improved = true;
      }
    }
    for (int sgn : {-1, 1}) {
      ShowerModelParams cand = best;
      cand.v2 = std::clamp(best.v2 + static_cast<double>(sgn) * step_v2, 1e-3, 2.0);
      const double c = eval_chi2(cand);
      if (c < best_chi2) {
        best = cand;
        best_chi2 = c;
        improved = true;
      }
    }
    for (int sgn : {-1, 1}) {
      ShowerModelParams cand = best;
      cand.energy = std::clamp(best.energy + static_cast<double>(sgn) * step_e, 0.2 * ee, 2.5 * ee);
      const double c = eval_chi2(cand);
      if (c < best_chi2) {
        best = cand;
        best_chi2 = c;
        improved = true;
      }
    }
    if (!improved) {
      step_v1 *= 0.5;
      step_v2 *= 0.5;
      step_e *= 0.5;
    }
  }

  return ShowerFitResult{.chi2 = best_chi2, .ndof = 8};
}

double compute_sigma_z_cm(float momentum_gev, float z_hpc_cm) noexcept {
  // 94D ring-dependent sigma_z parameterization (ELPDZS lines 12189-12226)
  constexpr std::array<double, 3> kA = {1.492, 2.048, 1.474};
  constexpr std::array<double, 3> kB = {-0.804, -0.865, -0.410};
  constexpr std::array<double, 3> kC = {0.251, 0.444, 0.431};

  const double za = std::fabs(static_cast<double>(z_hpc_cm));
  std::size_t ring = 0;
  if (za > 180.0) {
    ring = 2;
  } else if (za > 90.0) {
    ring = 1;
  }
  const double p = std::max(static_cast<double>(momentum_gev), 1.0);
  return 0.8 * std::sqrt(std::pow(kA[ring] * std::pow(p, kB[ring]), 2.0) +
                         kC[ring] * kC[ring]);
}

double compute_e_over_p_prob(double e_over_p, double momentum_gev) noexcept {
  if (!(e_over_p > 0.0) || e_over_p > 4.0) {
    return 0.0;
  }
  // DELPHI HPC resolution ~ 32%/sqrt(E) (+) 4.3%, with radiative low-E/p allowance
  const double e_est = std::max(momentum_gev, 1.0);
  double sigma = std::sqrt((0.32 * 0.32) / e_est + 0.043 * 0.043 + 0.015 * 0.015);
  if (e_over_p < 1.0) {
    // Allow modest bremsstrahlung tail on the low side, as in ELEVPN
    sigma *= 1.4;
  }
  const double chi2 = std::pow((e_over_p - 1.0) / sigma, 2.0);
  return detail::chi2_prob_1dof(std::min(chi2, 150.0));
}

}  // namespace

std::array<float, 9> expected_hpc_em_profile(float energy_gev,
                                             float theta_rad,
                                             float z_entry_cm) noexcept {
  std::array<float, 9> out{};
  if (!(energy_gev > 0.0f)) {
    return out;
  }
  const double ee = std::max(static_cast<double>(energy_gev), 1.0);
  const double log_ee = std::log(ee);
  const double z_abs = std::fabs(static_cast<double>(z_entry_cm));
  const ShowerModelParams params{
      .v1 = 1.0 / (4.234 + 0.7879 * log_ee),
      .v2 = (1.0 + 0.203e-3 * z_abs) / (2.8453 + 0.3042 * log_ee),
      .energy = static_cast<double>(energy_gev),
  };
  std::array<double, 9> raw{};
  compute_expected_9layers(static_cast<double>(theta_rad),
                           static_cast<double>(z_entry_cm),
                           params, raw);
  const double sum = std::accumulate(raw.begin(), raw.end(), 0.0);
  const double scale = (sum > 0.0) ? (static_cast<double>(energy_gev) / sum) : 1.0;
  for (std::size_t i = 0; i < 9; ++i) {
    out[i] = static_cast<float>(raw[i] * scale);
  }
  return out;
}

ElectronIdResult identify_electron(const ElectronIdInput& input) noexcept {
  ElectronIdResult res{};

  // ELECID requires charged track with p >= ELPMIN (1.0 GeV/c)
  if (input.charge == 0 || !(input.momentum_gev >= 1.0f)) {
    return res;
  }

  const double p = static_cast<double>(input.momentum_gev);
  const double ecal = std::max(static_cast<double>(input.ecal_energy_gev), 0.0);
  const double eop = ecal / p;
  res.e_over_p = static_cast<float>(eop);

  // Track-to-HPC spatial match (ELTAGC lines 10881-10885)
  const double sigma_z = compute_sigma_z_cm(input.momentum_gev, input.hpc_z_cm);
  bool match_ok = true;
  if (ecal > 0.5) {
    if (eop > 4.0 ||
        std::fabs(static_cast<double>(input.delta_z_cm)) > 4.0 * sigma_z ||
        std::fabs(static_cast<double>(input.delta_phi_pos_rad)) > 0.03 ||
        std::fabs(static_cast<double>(input.delta_phi_dir_rad)) > 0.5) {
      match_ok = false;
    }
  }
  res.passes_track_match = match_ok;
  if (!match_ok) {
    return res;
  }

  // Track-matching probabilities
  const double chi2_dz = std::pow(static_cast<double>(input.delta_z_cm) / sigma_z, 2.0);
  res.prob_dz = static_cast<float>(detail::chi2_prob_1dof(std::min(chi2_dz, 150.0)));

  const double sigma_dphi = 0.004 + 0.025 / p;
  const double chi2_dphi =
      std::pow(static_cast<double>(input.delta_phi_dir_rad) / sigma_dphi, 2.0);
  res.prob_dphi = static_cast<float>(detail::chi2_prob_1dof(std::min(chi2_dphi, 150.0)));

  // E/p probability
  if (ecal > 0.5) {
    res.prob_e_over_p = static_cast<float>(compute_e_over_p_prob(eop, p));
  }

  // Longitudinal HPC shower profile chi2
  const double layer_sum = std::accumulate(
      input.hpc_layer_energy_gev.begin(), input.hpc_layer_energy_gev.end(), 0.0);
  bool has_layers = false;
  if (layer_sum > 0.5) {
    has_layers = true;
    const float e_for_fit =
        (input.ecal_energy_gev > 0.5f) ? input.ecal_energy_gev : static_cast<float>(layer_sum);
    const auto fit = fit_hpc_longitudinal_profile(
        input.hpc_layer_energy_gev, e_for_fit, input.theta_rad,
        input.hpc_z_cm, input.is_crack);
    res.hpc_profile_chi2 = static_cast<float>(fit.chi2);
    res.hpc_ndof = fit.ndof;
    if (fit.chi2 >= 0.0 && fit.chi2 <= 150.0 && fit.ndof > 0) {
      res.prob_shower_fit =
          static_cast<float>(detail::chi2_prob(fit.chi2, fit.ndof));
    } else {
      res.prob_shower_fit = 0.0f;
    }
  } else if (ecal > 0.5) {
    // No layer breakdown supplied: keep neutral shower-shape fallback driven by E/p & HCAL
    res.prob_shower_fit = (input.hcal_energy_gev < 1.5f) ? 0.5f : 0.01f;
  }

  // TPC dE/dx pull for electron and pion (ELGDEDX)
  bool has_dedx = false;
  if (input.tpc_dedx > 0.0f && input.tpc_num_wires > 0) {
    has_dedx = true;
    const auto exp_e = bethe_bloch_dedx(input.momentum_gev, rich::Species::Electron);
    const auto exp_pi = bethe_bloch_dedx(input.momentum_gev, rich::Species::Pion);
    const float err_e = expected_dedx_error(exp_e.dedx_mip, input.tpc_num_wires, input.tpc_gap_cm);
    const float err_pi = expected_dedx_error(exp_pi.dedx_mip, input.tpc_num_wires, input.tpc_gap_cm);
    const float pull_e = gaussianize_dedx_pull((input.tpc_dedx - exp_e.dedx_mip) / err_e);
    const float pull_pi = gaussianize_dedx_pull((input.tpc_dedx - exp_pi.dedx_mip) / err_pi);
    res.dedx_pull_electron = pull_e;
    res.dedx_pull_pion = pull_pi;
    res.prob_dedx_electron = static_cast<float>(
        detail::chi2_prob_1dof(std::min(static_cast<double>(pull_e) * pull_e, 150.0)));
    res.prob_dedx_pion = static_cast<float>(
        detail::chi2_prob_1dof(std::min(static_cast<double>(pull_pi) * pull_pi, 150.0)));
  }

  // RICH electron probability
  bool has_rich = false;
  if (input.rich_electron_likelihood >= 0.0f) {
    has_rich = true;
    if (input.rich_pion_likelihood >= 0.0f) {
      const float sum = input.rich_electron_likelihood + input.rich_pion_likelihood;
      res.prob_rich_electron =
          (sum > 0.0f) ? (input.rich_electron_likelihood / sum) : 0.5f;
    } else {
      res.prob_rich_electron = std::clamp(input.rich_electron_likelihood, 0.0f, 1.0f);
    }
  }

  // Combined HPC & overall probabilities (ELPRTOT lines 10977-10989)
  if (res.prob_shower_fit > 1e-10f && res.prob_e_over_p > 0.05f) {
    const double chi2_hpc =
        detail::inverse_chi2_1dof(res.prob_e_over_p) +
        detail::inverse_chi2_1dof(res.prob_shower_fit) +
        detail::inverse_chi2_1dof(res.prob_dz) +
        detail::inverse_chi2_1dof(res.prob_dphi);
    res.prob_hpc = static_cast<float>(detail::chi2_prob(chi2_hpc, 4));

    const double dedx_pr = has_dedx
        ? static_cast<double>(res.prob_dedx_electron)
        : (0.43 / (1.0 + std::exp(-(std::min(80.0, p) - 11.0))));
    double chi2_all = chi2_hpc + detail::inverse_chi2_1dof(dedx_pr);
    int ndof_all = 5;
    if (has_rich) {
      chi2_all += detail::inverse_chi2_1dof(
          std::clamp(static_cast<double>(res.prob_rich_electron), 1e-6, 1.0));
      ++ndof_all;
    }
    res.prob_combined = static_cast<float>(detail::chi2_prob(chi2_all, ndof_all));
  }

  // HCAL hadronic leakage check
  const bool hcal_clean =
      (input.hcal_energy_gev < 2.0f) || (input.hcal_energy_gev < 0.35f * input.momentum_gev);

  // VeryLoose (ELVLOS lines 14717-14720)
  const bool very_loose =
      (eop > 0.3) || (has_dedx && input.tpc_dedx > 1.3f && res.prob_dedx_electron > 0.05f) ||
      (has_rich && res.prob_rich_electron > 0.55f);
  if (!very_loose) {
    return res;
  }
  res.tag = ElectronTag::VeryLoose;

  // Loose (ELTAGC lines 10895-10915)
  const bool eop_electron =
      (eop > std::min(0.90, std::max(0.5, 1.0 - 1.5 * 0.4 / std::sqrt(p)))) && (eop < 3.0);
  const bool dedx_or_fallback = has_dedx ? (res.prob_dedx_electron > 0.01f) : (eop > 0.65);
  const bool loose =
      ((dedx_or_fallback && res.prob_hpc > 1e-6f) ||
       (has_dedx && res.prob_dedx_electron > 0.10f && eop_electron) ||
       (res.prob_hpc > 0.01f && has_rich && res.prob_rich_electron > 0.65f) ||
       (has_dedx && res.prob_dedx_electron > 0.01f && has_rich && res.prob_rich_electron > 0.65f &&
        eop > 0.20));
  if (!loose || !hcal_clean) {
    return res;
  }
  res.tag = ElectronTag::Loose;

  // Standard (ELTAGC lines 10917-10941)
  if (has_layers && res.prob_shower_fit < 0.001f) {
    return res;
  }
  const bool std_dedx_ok = has_dedx ? (res.prob_dedx_electron > 0.002f) : (eop > 0.75);
  const bool standard =
      ((std_dedx_ok && res.prob_hpc > 0.01f) ||
       (res.prob_hpc > 0.10f && has_rich && res.prob_rich_electron > 0.75f) ||
       (has_dedx && res.prob_dedx_electron > 0.05f && has_rich && res.prob_rich_electron > 0.75f));
  if (!standard) {
    return res;
  }
  res.tag = ElectronTag::Standard;

  // Tight (ELTAGC lines 10943-10958)
  if (has_layers && res.prob_shower_fit < 0.005f) {
    return res;
  }
  if (input.is_conversion_veto) {
    return res;
  }
  const bool tight =
      (res.prob_combined > 0.15f) ||
      (res.prob_hpc > 0.15f && has_rich && res.prob_rich_electron > 0.85f) ||
      (has_dedx && res.prob_dedx_electron > 0.10f && has_rich && res.prob_rich_electron > 0.85f);
  if (tight) {
    res.tag = ElectronTag::Tight;
  }

  return res;
}

}  // namespace delphipp::reco::pid
