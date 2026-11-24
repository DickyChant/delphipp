// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Barrel & Forward RICH (BRICH / FRICH) particle ID
// (simana_v94d/src/car/richall.car: liquid C6F14 & gas C5F12/C4F10 Cherenkov
// angle reconstruction and 5-hypothesis mass likelihood for e, mu, pi, K, p).

#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::rich {

inline constexpr std::int32_t kDetectorIdBrich = 7;   // Barrel RICH
inline constexpr std::int32_t kDetectorIdFrich = 27;  // Forward RICH

inline constexpr float kRefractiveIndexLiquid = 1.273f;    // C6F14
inline constexpr float kRefractiveIndexGasBarrel = 1.00194f; // C5F12
inline constexpr float kRefractiveIndexGasFwd = 1.00153f;    // C4F10

enum class Species : std::size_t { Electron = 0, Muon = 1, Pion = 2, Kaon = 3, Proton = 4 };

// Particle rest masses [GeV/c^2]
inline constexpr std::array<float, 5> kSpeciesMassGeV = {
    0.000511f, 0.105658f, 0.139570f, 0.493677f, 0.938272f,
};

struct RichRingMeasurement {
  float momentum_gev = 0.0f;
  float refractive_index = kRefractiveIndexLiquid;
  float measured_theta_rad = 0.0f;
  float sigma_theta_rad = 0.005f;
  std::int32_t num_photons = 0;
};

struct RichLikelihoodResult {
  std::array<float, 5> expected_theta_rad{};  // -1 if below Cherenkov threshold
  std::array<float, 5> chi2{};
  Species best_species = Species::Pion;
};

// Expected Cherenkov angle theta_c = arccos(1 / (n * beta)), or nullopt below threshold.
[[nodiscard]] std::optional<float> expected_cherenkov_angle(
    float momentum_gev, float mass_gev, float refractive_index) noexcept;

// Evaluate 5-species Cherenkov hypothesis chi2 and select best mass hypothesis.
[[nodiscard]] RichLikelihoodResult evaluate_hypotheses(
    const RichRingMeasurement& ring) noexcept;

TanagraId write_ter_bank(TanagraStore& store,
                         const RichRingMeasurement& ring,
                         const RichLikelihoodResult& res,
                         std::int32_t detector_id = kDetectorIdBrich);

}  // namespace delphipp::reco::rich
