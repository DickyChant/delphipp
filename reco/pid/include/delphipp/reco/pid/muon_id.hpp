// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI muon identification
// (simana_v94d/src/car/muflag.car: MUFLAG / RTAG 94C/94D MUB/MUF/MUS hit chi2
// matching, plus simana_v94d/src/car/pxhaid.car: PXHAID / PXHAVE HCAL MIP
// energy veto) -> None / VeryLoose / Loose / Standard / Tight.

#pragma once

#include <array>
#include <cstdint>
#include <span>

namespace delphipp::reco::pid {

enum class MuonTag : std::uint8_t {
  None = 0,
  VeryLoose = 1,
  Loose = 2,
  Standard = 3,
  Tight = 4,
};

enum class MuonDetector : std::uint8_t {
  None = 0,
  Mub = 14,  // Barrel Muon Chambers (MUB, layers 1..7; outer > 3)
  Mus = 17,  // Surround Muon Chambers (SMC / MUS)
  Muf = 30,  // Forward Muon Chambers (MUF, layers 1..4; outer > 2)
};

struct MuonIdInput {
  float momentum_gev = 0.0f;
  float theta_rad = 1.5707963f;
  MuonDetector detector = MuonDetector::Mub;
  std::int32_t num_layers = 0;         // Associated muon chamber layers (MUTA(2))
  std::uint32_t hit_pattern = 0;       // Bitmask: bit (layer-1), bit 7 = overlap (MUTA(3))
  float chi2_total_pdf = 999.0f;       // Total fit chi2 / ndof (CHIPDF)
  float chi2_extrap_pdf = 999.0f;      // Track extrapolation chi2 / ndof (CHIEXT)
  // HCAL 4 longitudinal layer energies [GeV] (EH(1..4) in pxhaid.car)
  std::array<float, 4> hcal_layer_energy_gev{};
  // Fallback total HCAL energy [GeV] when layer breakdown is zero
  float hcal_total_energy_gev = 0.0f;
  bool apply_hcal_mip_veto = true;
};

struct MuonIdResult {
  MuonTag tag = MuonTag::None;
  bool has_outer_hit = false;
  float hcal_ecut_gev = 0.0f;
  float hcal_veto_rv = 0.0f;           // PXHAVE RV (< 1.0 = MIP consistent)
  bool passes_hcal_mip_veto = true;
  bool passes_muflag_veryloose = false;
  bool passes_muflag_loose = false;
  bool passes_muflag_standard = false;
  bool passes_muflag_tight = false;
};

// Compute (chi2_total_pdf, chi2_extrap_pdf) from per-layer chamber hit residuals
// and track-extrapolation chi2 (MUFLAG lines 1149-1178).
struct MuonChi2Summary {
  float chi2_total_pdf = 999.0f;
  float chi2_extrap_pdf = 999.0f;
};

[[nodiscard]] MuonChi2Summary compute_muon_hit_chi2(
    std::span<const float> hit_residuals_cm,
    std::span<const float> hit_errors_cm,
    float extrap_chi2,
    std::int32_t extrap_ndof) noexcept;

// Port of PXHAVE (pxhaid.car lines 224-320): compute HCAL MIP cut ECUT and veto RV.
[[nodiscard]] float compute_hcal_mip_veto(
    float momentum_gev,
    float theta_rad,
    const std::array<float, 4>& hcal_layers_gev,
    float hcal_total_energy_gev,
    float* ecut_gev_out = nullptr) noexcept;

// Evaluate MUFLAG + PXHAID muon identification.
[[nodiscard]] MuonIdResult identify_muon(const MuonIdInput& input) noexcept;

}  // namespace delphipp::reco::pid
