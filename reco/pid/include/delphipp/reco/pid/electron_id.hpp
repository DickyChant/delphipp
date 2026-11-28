// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI electron identification
// (simana_v94d/src/car/elepha32.car: ELECID / ELPRON / ELSHFI / ELPRTOT / ELTAGC):
// multi-variable e/pi discriminant combining E/p, HPC longitudinal shower-profile
// chi2, track-to-shower spatial matching (dz, dphi), TPC dE/dx pull, and RICH
// electron likelihood into None / VeryLoose / Loose / Standard / Tight tags.

#pragma once

#include <array>
#include <cstdint>

namespace delphipp::reco::pid {

enum class ElectronTag : std::uint8_t {
  None = 0,
  VeryLoose = 1,
  Loose = 2,
  Standard = 3,
  Tight = 4,
};

struct ElectronIdInput {
  float momentum_gev = 0.0f;
  std::int32_t charge = +1;
  float theta_rad = 1.5707963f;
  float ecal_energy_gev = 0.0f;
  float hcal_energy_gev = 0.0f;
  // 9 longitudinal HPC sampling-layer energies [GeV] (matching reco::hpc::kNumLayers).
  // Leave all zeros if longitudinal profile is not available.
  std::array<float, 9> hpc_layer_energy_gev{};
  float hpc_z_cm = 0.0f;            // Entry/shower z at HPC [cm]
  float delta_z_cm = 0.0f;          // z(extrap) - z(HPC) [cm]
  float delta_phi_pos_rad = 0.0f;   // phi(extrap) - phi(HPC) position [rad]
  float delta_phi_dir_rad = 0.0f;   // phi(extrap) - phi(HPC) direction [rad]
  float tpc_dedx = 0.0f;            // Truncated-mean dE/dx [MIP]; <= 0 if unavailable
  std::int32_t tpc_num_wires = 0;   // Number of TPC wires used
  float tpc_gap_cm = 0.4f;
  // RICH electron/pion likelihoods in [0, 1]; negative if RICH unavailable
  float rich_electron_likelihood = -1.0f;
  float rich_pion_likelihood = -1.0f;
  bool is_crack = false;            // HPC phi/theta crack region
  bool is_conversion_veto = false;  // Photon conversion veto (ELMACO IVETO > 0)
};

struct ElectronIdResult {
  ElectronTag tag = ElectronTag::None;
  float e_over_p = 0.0f;
  float hpc_profile_chi2 = -1.0f;
  std::int32_t hpc_ndof = 0;
  float dedx_pull_electron = 0.0f;
  float dedx_pull_pion = 0.0f;
  float prob_e_over_p = 0.0f;
  float prob_shower_fit = 0.0f;
  float prob_dz = 0.0f;
  float prob_dphi = 0.0f;
  float prob_dedx_electron = 0.0f;
  float prob_dedx_pion = 1.0f;
  float prob_rich_electron = -1.0f;
  float prob_hpc = 0.0f;       // Combined HPC probability (ELPRTOT PRHPC)
  float prob_combined = 0.0f;  // Combined e/pi discriminant probability (PROBEL)
  bool passes_track_match = false;
};

// Expected 9-layer HPC electromagnetic shower profile for an electron of given
// energy [GeV], polar angle [rad], and HPC entry z [cm] (ELSFMO/ELSIGMA).
[[nodiscard]] std::array<float, 9> expected_hpc_em_profile(
    float energy_gev, float theta_rad, float z_entry_cm = 0.0f) noexcept;

// Evaluate ELEPHA multi-variable electron identification.
[[nodiscard]] ElectronIdResult identify_electron(
    const ElectronIdInput& input) noexcept;

}  // namespace delphipp::reco::pid
