// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Forward Chambers A & B (FCA, FCB) reconstruction
// (simana_v94d/src/car/fcaana.car, fcbana.car): multi-orientation wire planes
// (u = x*cos(alpha) + y*sin(alpha)) solved via least squares into (x, y, z)
// track elements.

#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::fwd {

inline constexpr std::int32_t kDetectorIdFca = 20;  // FCA (|z| ~ 160 cm)
inline constexpr std::int32_t kDetectorIdFcb = 25;  // FCB (|z| ~ 275 cm)

struct ForwardWireHit {
  std::int32_t endcap = 1;  // 1 (+z) or 2 (-z)
  float z_cm = 160.0f;
  float alpha_rad = 0.0f;   // Measurement direction angle (0, +60, -60 deg)
  float u_cm = 0.0f;        // Measured coordinate u = x*cos(alpha) + y*sin(alpha)
  float sigma_u_cm = 0.03f;
};

struct ForwardSegment {
  std::int32_t endcap = 1;
  float x_cm = 0.0f;
  float y_cm = 0.0f;
  float z_cm = 0.0f;
  float chi2 = 0.0f;
  std::int32_t num_planes = 0;
};

// Solve (x, y) from >= 2 non-collinear wire-plane measurements via weighted least squares.
[[nodiscard]] std::optional<ForwardSegment> fit_segment(
    std::span<const ForwardWireHit> hits) noexcept;

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const ForwardSegment> segments,
    std::int32_t detector_id = kDetectorIdFca);

}  // namespace delphipp::reco::fwd
