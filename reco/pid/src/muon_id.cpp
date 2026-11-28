// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/pid/muon_id.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace delphipp::reco::pid {

namespace {

[[nodiscard]] bool check_outer_hit(MuonDetector det,
                                   std::uint32_t hit_pattern,
                                   std::int32_t num_layers) noexcept {
  // Bit 7 (1-based bit 8 in RTAG line 1332) marks overlap chambers
  if ((hit_pattern & (1u << 7)) != 0u) {
    return true;
  }
  switch (det) {
    case MuonDetector::Mub:
      // MUB outer/peripheral layers: 1-based J > 3 -> 0-based bits 3..6
      return (hit_pattern & 0b0111'1000u) != 0u;
    case MuonDetector::Muf:
      // MUF outer layers: 1-based J > 2 -> 0-based bits 2..3
      return (hit_pattern & 0b0000'1100u) != 0u;
    case MuonDetector::Mus:
      // SMC/MUS: any layer hit J >= 1
      return ((hit_pattern & 0b0111'1111u) != 0u) || (num_layers > 0);
    case MuonDetector::None:
      return false;
  }
  return false;
}

}  // namespace

MuonChi2Summary compute_muon_hit_chi2(std::span<const float> hit_residuals_cm,
                                      std::span<const float> hit_errors_cm,
                                      float extrap_chi2,
                                      std::int32_t extrap_ndof) noexcept {
  MuonChi2Summary out{};
  const std::size_t n = std::min(hit_residuals_cm.size(), hit_errors_cm.size());
  if (n == 0 || extrap_ndof <= 0) {
    return out;
  }

  double hit_chi2 = 0.0;
  std::size_t valid_hits = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double err = static_cast<double>(hit_errors_cm[i]);
    if (err > 1e-6) {
      const double pull = static_cast<double>(hit_residuals_cm[i]) / err;
      hit_chi2 += pull * pull;
      ++valid_hits;
    }
  }
  if (valid_hits == 0) {
    return out;
  }

  const double tot_ndof = static_cast<double>(valid_hits) + static_cast<double>(extrap_ndof);
  out.chi2_extrap_pdf =
      static_cast<float>(std::max(0.0, static_cast<double>(extrap_chi2)) /
                         static_cast<double>(extrap_ndof));
  out.chi2_total_pdf =
      static_cast<float>((hit_chi2 + std::max(0.0, static_cast<double>(extrap_chi2))) /
                         tot_ndof);
  return out;
}

float compute_hcal_mip_veto(float momentum_gev,
                            float theta_rad,
                            const std::array<float, 4>& hcal_layers_gev,
                            float hcal_total_energy_gev,
                            float* ecut_gev_out) noexcept {
  constexpr double kEca = 1.25;  // pxhaid.car line 48
  const double th_deg =
      std::clamp(static_cast<double>(theta_rad) * 180.0 / std::numbers::pi, 1.0, 179.0);
  const double sin_th = std::max(std::fabs(std::sin(static_cast<double>(theta_rad))), 0.1);
  const double x0 = std::sin(50.0 * std::numbers::pi / 180.0);

  double ecut = 0.0;
  if (th_deg > 50.0 && th_deg < 130.0) {
    ecut = kEca / (sin_th * sin_th);
  } else {
    ecut = kEca * (1.0 + sin_th) / ((1.0 + x0) * x0 * x0);
  }
  const double p_clamped = std::max(static_cast<double>(momentum_gev), 1.5);
  const double pcor = 1.0 + 0.25 * std::log(p_clamped);
  ecut *= pcor;

  if (ecut_gev_out != nullptr) {
    *ecut_gev_out = static_cast<float>(ecut);
  }

  double layer_sum = 0.0;
  for (float e : hcal_layers_gev) {
    layer_sum += std::max(0.0, static_cast<double>(e));
  }

  double rv = 0.0;
  if (layer_sum > 0.0) {
    const double half_cut = 0.5 * ecut;
    for (float e : hcal_layers_gev) {
      const double val = static_cast<double>(e);
      if (val > half_cut) {
        rv += (val / half_cut) - 1.0;
      }
    }
  } else if (hcal_total_energy_gev > 0.0f) {
    // Fallback when only total HCAL energy is available:
    // MIP deposits ~1-2 GeV total across 4 layers (well below 2 * ecut)
    const double eff_thresh = 2.0 * ecut;
    if (static_cast<double>(hcal_total_energy_gev) > eff_thresh) {
      rv = (static_cast<double>(hcal_total_energy_gev) / eff_thresh) - 1.0;
    }
  }

  return static_cast<float>(rv);
}

MuonIdResult identify_muon(const MuonIdInput& input) noexcept {
  MuonIdResult res{};

  float ecut = 0.0f;
  res.hcal_veto_rv = compute_hcal_mip_veto(
      input.momentum_gev, input.theta_rad,
      input.hcal_layer_energy_gev, input.hcal_total_energy_gev, &ecut);
  res.hcal_ecut_gev = ecut;
  res.passes_hcal_mip_veto = (res.hcal_veto_rv < 1.0f);

  if (input.detector == MuonDetector::None || input.num_layers <= 0 ||
      !(input.momentum_gev > 0.5f)) {
    return res;
  }

  res.has_outer_hit = check_outer_hit(input.detector, input.hit_pattern, input.num_layers);

  // JTAG = 4: VeryLoose (RTAG line 1675: at least 1 associated layer hit)
  res.passes_muflag_veryloose = (input.num_layers > 0);
  if (res.passes_muflag_veryloose) {
    res.tag = MuonTag::VeryLoose;
  }

  // Require >= 1.5 GeV/c for penetration cuts (pxhaid line 88)
  if (input.momentum_gev < 1.5f ||
      !std::isfinite(input.chi2_total_pdf) ||
      !std::isfinite(input.chi2_extrap_pdf) ||
      input.chi2_total_pdf < 0.0f || input.chi2_extrap_pdf < 0.0f) {
    return res;
  }

  // 94C/94D MUFLAG chi2 thresholds (muflag.car RTAG lines 1356-1667)
  float loose_tot = 0.0f;
  float loose_ext = 0.0f;
  float std_tot = 0.0f;
  float std_ext = 0.0f;
  float tight_tot = 0.0f;
  float tight_ext = 0.0f;

  switch (input.detector) {
    case MuonDetector::Mub:
      loose_tot = 6.0f;  loose_ext = 6.0f;
      std_tot = 3.5f;    std_ext = 5.25f;
      tight_tot = 2.0f;  tight_ext = 3.0f;
      break;
    case MuonDetector::Muf:
      loose_tot = 7.0f;  loose_ext = 7.0f;
      std_tot = 4.0f;    std_ext = 3.0f;
      tight_tot = 2.0f;  tight_ext = 1.5f;
      break;
    case MuonDetector::Mus:
      loose_tot = 6.0f;  loose_ext = 5.0f;
      std_tot = 4.0f;    std_ext = 5.0f;
      tight_tot = 2.0f;  tight_ext = 5.0f;
      break;
    case MuonDetector::None:
      return res;
  }

  res.passes_muflag_loose =
      (input.chi2_total_pdf <= loose_tot && input.chi2_extrap_pdf <= loose_ext);
  res.passes_muflag_standard =
      res.has_outer_hit &&
      (input.chi2_total_pdf <= std_tot && input.chi2_extrap_pdf <= std_ext);
  res.passes_muflag_tight =
      res.has_outer_hit &&
      (input.chi2_total_pdf <= tight_tot && input.chi2_extrap_pdf <= tight_ext);

  const bool veto_ok = !input.apply_hcal_mip_veto || res.passes_hcal_mip_veto;

  if (res.passes_muflag_loose && veto_ok) {
    res.tag = MuonTag::Loose;
  }
  if (res.passes_muflag_standard && veto_ok) {
    res.tag = MuonTag::Standard;
  }
  if (res.passes_muflag_tight && veto_ok) {
    res.tag = MuonTag::Tight;
  }

  return res;
}

}  // namespace delphipp::reco::pid
