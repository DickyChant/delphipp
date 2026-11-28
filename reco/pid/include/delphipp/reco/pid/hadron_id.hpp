// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI hadron identification:
//  - simana_v94d/src/car/dedxid.car (GETDNEW Bethe-Bloch ionization curve,
//    BBDXER wire/gap error model, BBPULL Gaussianization, BBDXID likelihoods)
//  - simana_v94d/src/car/pxrich.car & richall.car (5-species RICH Cherenkov
//    combined likelihood for e, mu, pi, K, p)
//  - simana_v94d/src/car/pxhaid.car (hadron vs MIP context)

#pragma once

#include <array>
#include <cstdint>
#include <optional>

#include "delphipp/reco/rich.hpp"

namespace delphipp::reco::pid {

struct DedxExpectation {
  float dedx_mip = -1.0f;  // Expected dE/dx in minimum-ionizing units (~1.0 at MIP)
  float ddedx_dp = -1.0f;  // Slope d(dE/dx)/dp [MIP / (GeV/c)]
};

struct HadronIdInput {
  float momentum_gev = 0.0f;
  // TPC 80% truncated-mean dE/dx measurement
  float tpc_dedx = -1.0f;         // Measured dE/dx [MIP]; <= 0 if unavailable
  std::int32_t tpc_num_wires = 0; // Number of TPC wires used (typically 40..192)
  float tpc_gap_cm = 0.4f;        // TPC wire gap normalization (~0.4 cm)
  bool tpc_edge_region = false;   // Sector edge flag (NMETH == 2 in BBDXER)
  // RICH ring measurements (liquid C6F14 and/or gas C5F12/C4F10)
  std::optional<rich::RichRingMeasurement> rich_liquid;
  std::optional<rich::RichRingMeasurement> rich_gas;
  // Optional precomputed 5-species RICH chi2 (order: e, mu, pi, K, p)
  std::optional<std::array<float, 5>> rich_chi2_override;
  bool gaussianize_dedx_pull = true;
  bool is_simulation = false;
};

struct HadronIdResult {
  std::array<float, 5> expected_dedx{};
  std::array<float, 5> dedx_error{};
  std::array<float, 5> dedx_pull{};
  std::array<float, 5> dedx_chi2{};
  std::array<float, 5> rich_chi2{};
  std::array<float, 5> total_chi2{};
  std::array<float, 5> likelihood{};  // Normalized 5-species posterior probabilities
  std::int32_t dedx_quality = -1;     // -1 = no info, 0 = poor .. 4 = perfect (BBDXER)
  rich::Species best_species = rich::Species::Pion;
  float prob_kaon_vs_pion = 0.5f;
  float prob_proton_vs_pion = 0.5f;
  float prob_proton_vs_kaon = 0.5f;
  std::uint8_t kaon_tag = 0;          // 0 = None, 1 = Loose, 2 = Standard, 3 = Tight
  std::uint8_t proton_tag = 0;        // 0 = None, 1 = Loose, 2 = Standard, 3 = Tight
};

// Exact port of DELPHI 94C/94D Bethe-Bloch curve GETDNEW (dedxid.car lines 4109-4190).
[[nodiscard]] DedxExpectation bethe_bloch_dedx(float momentum_gev,
                                               float mass_gev) noexcept;

[[nodiscard]] DedxExpectation bethe_bloch_dedx(float momentum_gev,
                                               rich::Species species) noexcept;

// Port of BBDXER (dedxid.car lines 1120-1186): expected dE/dx resolution & quality (-1..4).
[[nodiscard]] float expected_dedx_error(float expected_dedx_mip,
                                        std::int32_t num_wires,
                                        float gap_cm = 0.4f,
                                        bool edge_region = false,
                                        std::int32_t* quality_out = nullptr) noexcept;

// Port of BBPULL (dedxid.car lines 1512-1669): Landau-tail Gaussianization of dE/dx pull.
[[nodiscard]] float gaussianize_dedx_pull(float raw_pull,
                                          bool is_simulation = false) noexcept;

// Port of BBDXID (dedxid.car lines 310-429): pairwise TPC dE/dx separation probability.
[[nodiscard]] float dedx_pairwise_probability(float momentum_gev,
                                              float measured_dedx,
                                              std::int32_t num_wires,
                                              rich::Species signal_species,
                                              rich::Species background_species,
                                              float gap_cm = 0.4f,
                                              bool edge_region = false) noexcept;

// Combined TPC dE/dx + RICH (liquid & gas) 5-species identification (e, mu, pi, K, p).
[[nodiscard]] HadronIdResult identify_hadron(const HadronIdInput& input) noexcept;

}  // namespace delphipp::reco::pid
