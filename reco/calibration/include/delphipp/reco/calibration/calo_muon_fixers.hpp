// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Calorimeter, Muon Chamber, and RICH offline
// calibration fixers (simana_v94d/src/car/hacfix36.car, haccorxx.car,
// mucfix.car, ricfix.car).

#pragma once

#include <cstdint>
#include <optional>

namespace delphipp::reco::calibration {

// HCAL (HAC) tower quality status and per-tower gain equalization (hacfix36 / haccorxx).
enum class HcalTowerStatus : std::uint8_t {
  Good = 0,
  Dead = 1,
  Hot = 2
};

struct HcalTowerCalibration {
  std::int32_t itheta = 0;
  std::int32_t iphi = 0;
  std::int32_t ilayer = 1;
  HcalTowerStatus status = HcalTowerStatus::Good;
  float gain_factor = 1.0f;
};

// Environmental gas parameters matching ECALIB (hadana.car / hacfix36.car).
struct HcalFixConfig {
  std::int32_t run_number = 47014;  // Real data if > 10581
  float pressure_mbar = 975.0f;     // Active window (950, 1000) mbar
  float isobutane_pct = 30.0f;      // Active window (29, 31) %
  float temperature_c = 20.0f;      // Active window (16, 23) deg C
  float run_table_scale = 1.0f;     // HCV33/HCV34/HCV35 Z0 hadronic scale
  bool enforce_sensor_windows = true;
};

// Return HCV33 step-period run scale factor from hacfix36.car (default 1.0f).
[[nodiscard]] float hcal_run_table_scale(std::int32_t run_number) noexcept;

// Apply HCAL gas pressure/isobutane/temperature scaling, run-table scale, and
// dead/hot tower masking & gain equalization. Returns std::nullopt if masked/invalid.
[[nodiscard]] std::optional<float> correct_hcal_tower_energy(
    float raw_energy_gev,
    const HcalFixConfig& env = {},
    const HcalTowerCalibration& tower = {}) noexcept;

// Muon chamber subsystems covered by mucfix.car.
enum class MuonSubsystem : std::int32_t {
  Mub = 14,
  Mus = 17,
  Muf = 30
};

struct MuonChamberFix {
  MuonSubsystem subsystem = MuonSubsystem::Mub;
  std::int32_t chamber_id = 1;
  std::int32_t layer = 1;
  float rphi_offset_cm = 0.0f;      // DBRHIT / RFISFT residual shift [cm]
  float z_offset_cm = 0.0f;         // DBZHIT delay-line z residual shift [cm]
  float anode_error_scale = 1.0f;   // SCDBAN
  float delay_error_scale = 1.0f;   // SCDBCT
  float efficiency = 1.0f;          // Chamber efficiency in [0, 1]
  bool is_dead = false;
};

struct MuonChamberHit {
  MuonSubsystem subsystem = MuonSubsystem::Mub;
  std::int32_t chamber_id = 1;
  std::int32_t layer = 1;
  float rphi_cm = 0.0f;
  float z_cm = 0.0f;
  float sigma_rphi_cm = 0.2f;
  float sigma_z_cm = 1.0f;
  float weight = 1.0f;
};

// Apply muon chamber residual offsets, error scaling, and efficiency/dead-chamber filtering.
// If `efficiency_roll` in [0, 1) is supplied, hits with `*efficiency_roll >= fix.efficiency`
// are rejected (returning std::nullopt) to reproduce data efficiency in MC.
[[nodiscard]] std::optional<MuonChamberHit> correct_muon_chamber(
    const MuonChamberHit& hit,
    const MuonChamberFix& fix,
    std::optional<float> efficiency_roll = std::nullopt) noexcept;

// DELPHI RICH radiator media (ricfix.car / rich.hpp).
enum class RichRadiator : std::uint8_t {
  LiquidC6F14 = 0,
  GasBarrelC5F12 = 1,
  GasForwardC4F10 = 2
};

inline constexpr float kNominalIndexC6F14 = 1.273f;
inline constexpr float kNominalIndexC5F12 = 1.00194f;
inline constexpr float kNominalIndexC4F10 = 1.00153f;
inline constexpr float kDefaultRichRefPressureMbar = 1013.25f;
inline constexpr float kDefaultRichRefTempBarrelK = 313.15f;   // 40 C (BRICH)
inline constexpr float kDefaultRichRefTempForwardK = 298.15f;  // 25 C (FRICH)

// Gladstone-Dale P/T scaling of C6F14 (liquid), C5F12 (BRICH gas), and C4F10 (FRICH gas)
// refractive indices (ricfix.car).
// Pass `nominal_index <= 1.0f` or `ref_temperature_k <= 0.0f` to use the radiator defaults.
[[nodiscard]] float adjusted_rich_refractive_index(
    RichRadiator radiator,
    float pressure_mbar,
    float temperature_k,
    float nominal_index = 0.0f,
    float ref_pressure_mbar = kDefaultRichRefPressureMbar,
    float ref_temperature_k = 0.0f) noexcept;

// Adjust Cherenkov angle for updated refractive index: cos(theta_new) = (n_old / n_new) * cos(theta_old).
[[nodiscard]] float correct_cherenkov_angle(
    float measured_theta_rad,
    float n_old,
    float n_new) noexcept;

}  // namespace delphipp::reco::calibration
