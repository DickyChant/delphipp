// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 energy-flow, Durham (kT) jet clustering, and event-shape
// algorithms (replacing simana_v94d/src/car/eflowxx.car & vecsub72.car).

#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "delphipp/dst/short_dst.hpp"
#include "delphipp/geometry.hpp"

namespace delphipp::dst {

struct Jet {
  float px = 0.0f;
  float py = 0.0f;
  float pz = 0.0f;
  float energy = 0.0f;
  std::vector<std::size_t> constituent_indices;
};

struct EventShape {
  float thrust = 0.5f;
  Vec3 thrust_axis = {0.0f, 0.0f, 1.0f};
  float visible_mass_gev = 0.0f;
  float missing_pt_gev = 0.0f;
  float y_cut_3to2 = 0.0f;
};

// Cluster Short-DST particles with the e+e- Durham (kT) algorithm down to
// scale y_cut = 2 * min(E_i^2, E_j^2) * (1 - cos(theta_ij)) / E_vis^2.
[[nodiscard]] std::vector<Jet> cluster_durham_jets(
    std::span<const DstParticle> particles,
    float y_cut = 0.01f,
    std::size_t min_jets = 2);

// Compute event thrust, thrust axis, visible invariant mass, missing pT,
// and the Durham y_32 (3-to-2 jet) transition scale.
[[nodiscard]] EventShape compute_event_shape(
    std::span<const DstParticle> particles);

}  // namespace delphipp::dst
