// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 V0 & secondary decay vertex reconstruction
// (K_S^0 -> pi+ pi-, Lambda -> p pi-, Lambda_bar -> p_bar pi+, gamma -> e+ e-)
// replacing core of simana_v94d/src/car/recv022.car & rec2nd.car.

#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/geometry.hpp"
#include "delphipp/helix.hpp"

namespace delphipp::reco::tracking {

enum class V0Hypothesis : std::int32_t {
  PhotonConversion = 22,
  KShort = 310,
  Lambda = 3122,
  AntiLambda = -3122,
};

struct V0Candidate {
  std::size_t pos_idx = 0;
  std::size_t neg_idx = 0;
  Vec3 decay_xyz_cm{};
  float decay_length_cm = 0.0f;
  float dca_cm = 0.0f;
  float mass_kshort_gev = 0.0f;
  float mass_lambda_gev = 0.0f;
  float mass_antilambda_gev = 0.0f;
  float mass_ee_gev = 0.0f;
  V0Hypothesis best_hypothesis = V0Hypothesis::KShort;
};

[[nodiscard]] std::vector<V0Candidate> find_v0_candidates(
    std::span<const PerigeeTrack> tracks,
    std::span<const std::int32_t> charges,
    float bz_gev_per_cm,
    const Vec3& pv_xyz_cm = {0.0f, 0.0f, 0.0f},
    float max_dca_cm = 1.5f);

}  // namespace delphipp::reco::tracking
