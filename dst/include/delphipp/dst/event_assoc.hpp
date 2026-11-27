// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 MC truth-to-reconstructed particle association
// (replacing core of simana_v94d/src/car/evasso.car).

#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "delphipp/dst/short_dst.hpp"

namespace delphipp::dst {

struct TrackMatch {
  std::size_t reco_idx = 0;
  std::size_t truth_idx = 0;
  float delta_r = 0.0f;
  float rel_dp = 0.0f;
};

// Greedy 1-to-1 best angular + charge-consistent matching between reconstructed
// Short-DST particles and MC truth particles.
[[nodiscard]] std::vector<TrackMatch> associate_particles(
    std::span<const DstParticle> reco,
    std::span<const DstParticle> truth,
    float max_delta_r = 0.15f);

}  // namespace delphipp::dst
