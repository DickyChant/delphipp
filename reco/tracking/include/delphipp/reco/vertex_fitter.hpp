// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 Primary Vertex fitter with beamspot constraint (simana_v94d/src/car/vertex.car).

#pragma once

#include <array>
#include <cstdint>
#include <span>

#include "delphipp/geometry.hpp"
#include "delphipp/helix.hpp"

namespace delphipp::reco::tracking {

struct BeamspotConstraint {
  Vec3 xyz_cm = {0.0f, 0.0f, 0.0f};
  Vec3 sigma_cm = {0.015f, 0.001f, 0.70f};  // Typical LEP1 beamspot (150 um x 10 um x 7 mm)
};

struct VertexFitResult {
  Vec3 xyz_cm{};
  Vec3 sigma_cm{};
  float chi2 = 0.0f;
  std::int32_t ndf = 0;
  std::int32_t num_tracks = 0;
};

// Linearized weighted least-squares primary vertex fit from perigee tracks + beamspot.
[[nodiscard]] VertexFitResult fit_primary_vertex(
    std::span<const PerigeeTrack> tracks,
    const BeamspotConstraint& beamspot = {},
    float track_sigma_rphi_cm = 0.005f,
    float track_sigma_z_cm = 0.020f) noexcept;

}  // namespace delphipp::reco::tracking
