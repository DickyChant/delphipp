// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELSIM simulation core & fast GEANT3 stepping
// (simana_v94d/src/car/delsim36.car, geant315.car):
//   - SimParticle / SimEvent truth model
//   - Deterministic LEP1 Gaussian beamspot primary-vertex smearing
//   - RK4 Lorentz-force stepping across DELPHI tracker radii & endcap planes
//     with PDG Highland multiple scattering and Bethe-Bloch ionization loss.

#pragma once

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "delphipp/geometry.hpp"
#include "delphipp/magnetic_field.hpp"

namespace delphipp::sim {

struct SimParticle {
  std::int32_t pdg = 0;
  std::int32_t charge = 0;
  float px = 0.0f;
  float py = 0.0f;
  float pz = 0.0f;
  float energy = 0.0f;
  float mass = 0.0f;
  Vec3 vertex_cm{};
};

struct SimEvent {
  std::int32_t run = 0;
  std::int32_t event = 1;
  Vec3 primary_vertex_cm{};
  std::vector<SimParticle> particles;
};

struct BeamspotProfile {
  Vec3 nominal_cm{0.0f, 0.0f, 0.0f};
  float sigma_x_cm = 0.015f;
  float sigma_y_cm = 0.001f;
  float sigma_z_cm = 1.0f;
};

// Deterministic LEP1 Gaussian beamspot smearer backed by std::mt19937_64.
class PrimaryVertexSmearer {
 public:
  explicit PrimaryVertexSmearer(std::uint64_t seed = 5489ULL,
                                BeamspotProfile profile = {}) noexcept;
  PrimaryVertexSmearer(std::uint64_t seed,
                       float sigma_x_cm,
                       float sigma_y_cm,
                       float sigma_z_cm,
                       Vec3 nominal_cm = {}) noexcept;

  [[nodiscard]] Vec3 sample() noexcept;
  [[nodiscard]] const BeamspotProfile& profile() const noexcept { return profile_; }
  void reseed(std::uint64_t seed) noexcept { rng_.seed(seed); }

 private:
  BeamspotProfile profile_{};
  std::mt19937_64 rng_;
};

struct SteppedHit {
  std::string det_name;
  float r_cm = 0.0f;
  float phi_rad = 0.0f;
  float z_cm = 0.0f;
  float path_x0 = 0.0f;
  float dedx_mev = 0.0f;
};

// Step a particle through `bfield` across DELPHI tracker layers
// (VD 6.3/9.0/10.9 cm, ID 12..22 cm, STW 23.2..28 cm, TPC 35..110 cm,
// OD 197..206 cm, FCA |z|=160 cm, FCB |z|=275 cm) with Highland MS and
// Bethe-Bloch ionization energy loss.
[[nodiscard]] std::vector<SteppedHit> step_particle_in_bfield(
    const SimParticle& particle,
    const MagneticField& bfield,
    std::uint64_t seed = 1ULL);

}  // namespace delphipp::sim
