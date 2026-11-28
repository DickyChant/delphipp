// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI 3D lifetime b-tagging
// (simana_v94d/src/car/aabtagxx.car: AADIST / AASGNT / AAIMPC / AATRKP /
// bbPROB / AAEPRB / AAHEMI): signed 3D impact parameter significance S = d0/sigma
// wrt jet axis, calibrated resolution function P_T(S), and combined
// track / hemisphere / event b-tag probability P_E via Fisher log-likelihood.

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/geometry.hpp"
#include "delphipp/helix.hpp"
#include "delphipp/reco/vertex_fitter.hpp"

namespace delphipp::reco::pid {

struct BTagConfig {
  float min_momentum_gev = 0.2f;       // PAR9
  float max_impact_rphi_cm = 0.25f;    // PAR6 (2.5 mm)
  float max_impact_z_cm = 0.30f;       // PAR11 (3.0 mm)
  float max_dist_to_jet_sig = 10.0f;   // PAR12
  float min_sigma_cm = 1.0e-4f;        // 1 um guard
  float bz_gev_per_cm = 0.00360f;      // ~1.2 T DELPHI field (GeV/c / cm)
  bool use_3d_z = true;
};

struct BTagTrackInput {
  PerigeeTrack track{};
  float sigma_rphi_cm = 0.0035f;   // Perigee r-phi uncertainty [cm] (~35 um with VD)
  float sigma_z_cm = 0.0070f;      // Perigee z uncertainty [cm] (~70 um with VD z)
  std::int32_t n_vd_rphi = 3;      // Associated VD r-phi layers (1..3)
  std::int32_t n_vd_z = 2;         // Associated VD z layers (0..3)
  Vec3 jet_axis{};                 // Associated jet direction; falls back to thrust_axis if zero
};

struct BTagTrackResult {
  float ip_rphi_cm = 0.0f;         // Signed r-phi impact parameter [cm]
  float ip_z_cm = 0.0f;            // Signed z impact parameter [cm]
  float ip_3d_cm = 0.0f;           // Signed 3D impact parameter wrt PV [cm]
  float sigma_rphi_cm = 0.0f;
  float sigma_z_cm = 0.0f;
  float sigma_3d_cm = 0.0f;
  float significance_rphi = 0.0f;  // Signed S_rphi = sign * |d_rphi| / sigma_rphi
  float significance_z = 0.0f;     // Signed S_z = sign * |d_z| / sigma_z
  float significance_3d = 0.0f;    // Signed S_3D = sign * |d_3D| / sigma_3D
  float dist_to_jet_cm = 0.0f;     // 3D DCA distance from track to jet axis
  float prob_rphi = 1.0f;          // P_T(|S_rphi|)
  float prob_z = 1.0f;             // P_T(|S_z|)
  float prob_3d = 1.0f;            // P_T(|S_3D|)
  std::int32_t hemisphere = 1;     // 1 (along thrust) or 2 (opposite thrust)
  bool selected = false;
};

struct BTagEventResult {
  std::vector<BTagTrackResult> tracks;
  std::int32_t num_selected_tracks = 0;
  std::int32_t num_positive_tracks = 0;
  std::int32_t num_negative_tracks = 0;
  // Event b-tag probabilities (Fisher log-likelihood combination bbPROB):
  // Low p_event_btag (< 1e-3) indicates displaced heavy-flavour tracks (b-jet).
  float p_event_btag = 1.0f;       // PROBP (tracks with S >= 0)
  float p_event_negative = 1.0f;   // PROBN (tracks with S < 0, resolution control)
  float p_event_all = 1.0f;        // PROBS (all selected tracks)
  std::array<float, 2> p_hemisphere_btag{1.0f, 1.0f};      // PROBHP[0..1]
  std::array<float, 2> p_hemisphere_negative{1.0f, 1.0f};  // PROBHN[0..1]
};

// Calibrated DELPHI resolution probability P_T(S) for a prompt track to have
// impact parameter significance >= |S| (AATRKP parameterization).
// If custom_table_400 is non-empty (400 bins covering |S| in [0, 40)), performs
// exact AATRKP quadratic interpolation; otherwise uses the 94D multi-Gaussian+tail model.
[[nodiscard]] float track_resolution_probability(
    float abs_significance,
    std::int32_t n_vd_layers = 3,
    std::span<const float> custom_table_400 = {}) noexcept;

// Fisher log-likelihood combination of track p-values (exact port of bbPROB).
[[nodiscard]] float fisher_combine_probabilities(
    std::span<const float> track_probs) noexcept;

// Compute 3D signed impact parameter & significances for a single track wrt PV & jet axis.
[[nodiscard]] BTagTrackResult evaluate_track_btag(
    const BTagTrackInput& track_in,
    const reco::tracking::VertexFitResult& primary_vertex,
    const Vec3& thrust_axis,
    const BTagConfig& cfg = {}) noexcept;

// Run full AABTAG 3D lifetime b-tagging over all event tracks.
[[nodiscard]] BTagEventResult evaluate_event_btag(
    std::span<const BTagTrackInput> tracks,
    const reco::tracking::VertexFitResult& primary_vertex,
    const Vec3& thrust_axis,
    const BTagConfig& cfg = {}) noexcept;

}  // namespace delphipp::reco::pid
