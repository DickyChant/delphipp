// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/calibration/calo_muon_fixers.hpp"

#include <algorithm>
#include <cmath>

namespace delphipp::reco::calibration {

float hcal_run_table_scale(std::int32_t run_number) noexcept {
  // Matches hacfix36.car HCV33 period steps (1993 reprocessing).
  if (run_number >= 41600 && run_number <= 41632) {
    return (4.1609f + (5.2272f - 4.1609f) / 3.0f) / 4.1609f;
  }
  if (run_number >= 41633 && run_number <= 41647) {
    return (4.1609f + 2.0f * (5.2272f - 4.1609f) / 3.0f) / 4.1609f;
  }
  if (run_number >= 41650 && run_number <= 42033) {
    return 5.2272f / 4.1609f;
  }
  if (run_number == 42034) {
    return (5.2272f - (5.2272f - 4.5276f) / 3.0f) / 4.1609f;
  }
  if (run_number >= 42075 && run_number <= 43093) {
    return 4.5276f / 4.1609f;
  }
  return 1.0f;
}

std::optional<float> correct_hcal_tower_energy(
    float raw_energy_gev,
    const HcalFixConfig& env,
    const HcalTowerCalibration& tower) noexcept {
  if (tower.status == HcalTowerStatus::Dead || tower.status == HcalTowerStatus::Hot) {
    return std::nullopt;
  }
  if (!std::isfinite(raw_energy_gev) || raw_energy_gev <= 0.0f ||
      !std::isfinite(tower.gain_factor) || tower.gain_factor <= 0.0f ||
      !std::isfinite(env.run_table_scale) || env.run_table_scale <= 0.0f) {
    return std::nullopt;
  }

  double energy = static_cast<double>(raw_energy_gev) *
                  static_cast<double>(tower.gain_factor) *
                  static_cast<double>(env.run_table_scale);

  if (env.run_number > 10581) {
    if (std::isfinite(env.pressure_mbar) &&
        (!env.enforce_sensor_windows ||
         (env.pressure_mbar > 950.0f && env.pressure_mbar < 1000.0f))) {
      energy *= 1.0 + 0.008056 * (static_cast<double>(env.pressure_mbar) - 975.0);
    }
    if (std::isfinite(env.isobutane_pct) &&
        (!env.enforce_sensor_windows ||
         (env.isobutane_pct > 29.0f && env.isobutane_pct < 31.0f))) {
      energy *= 1.0 + 0.19 * (static_cast<double>(env.isobutane_pct) - 30.0);
    }
    if (std::isfinite(env.temperature_c) &&
        (!env.enforce_sensor_windows ||
         (env.temperature_c > 16.0f && env.temperature_c < 23.0f))) {
      energy *= 1.0 - 0.030 * (static_cast<double>(env.temperature_c) - 20.0);
    }
  }

  if (energy <= 0.0 || !std::isfinite(energy)) {
    return std::nullopt;
  }
  return static_cast<float>(energy);
}

std::optional<MuonChamberHit> correct_muon_chamber(
    const MuonChamberHit& hit,
    const MuonChamberFix& fix,
    std::optional<float> efficiency_roll) noexcept {
  if (fix.is_dead || !std::isfinite(fix.efficiency) || fix.efficiency <= 0.0f) {
    return std::nullopt;
  }
  const float eff = std::clamp(fix.efficiency, 0.0f, 1.0f);
  if (efficiency_roll.has_value()) {
    if (!std::isfinite(*efficiency_roll) || *efficiency_roll >= eff) {
      return std::nullopt;
    }
  }

  MuonChamberHit out = hit;
  out.rphi_cm += fix.rphi_offset_cm;
  out.z_cm += fix.z_offset_cm;
  if (std::isfinite(fix.anode_error_scale) && fix.anode_error_scale > 0.0f) {
    out.sigma_rphi_cm *= fix.anode_error_scale;
  }
  if (std::isfinite(fix.delay_error_scale) && fix.delay_error_scale > 0.0f) {
    out.sigma_z_cm *= fix.delay_error_scale;
  }
  out.weight *= eff;
  return out;
}

float adjusted_rich_refractive_index(RichRadiator radiator,
                                     float pressure_mbar,
                                     float temperature_k,
                                     float nominal_index,
                                     float ref_pressure_mbar,
                                     float ref_temperature_k) noexcept {
  double n0 = static_cast<double>(nominal_index);
  double t0 = static_cast<double>(ref_temperature_k);
  switch (radiator) {
    case RichRadiator::LiquidC6F14:
      if (n0 <= 1.0) n0 = static_cast<double>(kNominalIndexC6F14);
      if (t0 <= 0.0) t0 = static_cast<double>(kDefaultRichRefTempBarrelK);
      break;
    case RichRadiator::GasBarrelC5F12:
      if (n0 <= 1.0) n0 = static_cast<double>(kNominalIndexC5F12);
      if (t0 <= 0.0) t0 = static_cast<double>(kDefaultRichRefTempBarrelK);
      break;
    case RichRadiator::GasForwardC4F10:
      if (n0 <= 1.0) n0 = static_cast<double>(kNominalIndexC4F10);
      if (t0 <= 0.0) t0 = static_cast<double>(kDefaultRichRefTempForwardK);
      break;
  }

  const double p0 = (ref_pressure_mbar > 0.0f && std::isfinite(ref_pressure_mbar))
                        ? static_cast<double>(ref_pressure_mbar)
                        : static_cast<double>(kDefaultRichRefPressureMbar);

  if (!std::isfinite(pressure_mbar) || pressure_mbar <= 0.0f ||
      !std::isfinite(temperature_k) || temperature_k <= 0.0f) {
    return static_cast<float>(n0);
  }

  const double p = static_cast<double>(pressure_mbar);
  const double t = static_cast<double>(temperature_k);

  if (radiator == RichRadiator::LiquidC6F14) {
    // Liquid C6F14 Gladstone-Dale density scaling: thermal expansion beta ~ 1.6e-3 / K,
    // isothermal compressibility kappa_T ~ 1.5e-7 / mbar.
    constexpr double kBetaLiquidPerK = 1.6e-3;
    constexpr double kKappaLiquidPerMbar = 1.5e-7;
    const double rho_ratio =
        std::max(0.5, 1.0 - kBetaLiquidPerK * (t - t0) + kKappaLiquidPerMbar * (p - p0));
    return static_cast<float>(1.0 + (n0 - 1.0) * rho_ratio);
  }

  // Ideal-gas Gladstone-Dale scaling: (n - 1) / (n0 - 1) = (P / P0) * (T0 / T).
  const double density_ratio = (p / p0) * (t0 / t);
  return static_cast<float>(1.0 + (n0 - 1.0) * density_ratio);
}

float correct_cherenkov_angle(float measured_theta_rad,
                              float n_old,
                              float n_new) noexcept {
  if (!std::isfinite(measured_theta_rad) || !std::isfinite(n_old) ||
      !std::isfinite(n_new) || n_new <= 1.0f || n_old <= 1.0f) {
    return measured_theta_rad;
  }
  const double cos_new =
      (static_cast<double>(n_old) / static_cast<double>(n_new)) *
      std::cos(static_cast<double>(measured_theta_rad));
  return static_cast<float>(std::acos(std::clamp(cos_new, -1.0, 1.0)));
}

}  // namespace delphipp::reco::calibration
