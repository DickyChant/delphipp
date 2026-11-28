// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 DELSIM calorimeter, timing, muon, and RICH digitization
// (hpcsim, emfsim, hacsim, stcsim, lumsim, tofsim, hofsim, mubsim, mufsim,
// mussim, richsim) producing relocated ZEBRA raw banks ready for
// EventReconstructionPipeline.

#pragma once

#include <cstdint>
#include <span>

#include "delphipp/io/fz_zebra_parser.hpp"

namespace delphipp::sim {

struct SimShowerParticle {
  std::int32_t pdg = 0;
  std::int32_t charge = 0;
  float px = 0.0f;
  float py = 0.0f;
  float pz = 0.0f;
  float energy = 0.0f;
};

[[nodiscard]] io::ParsedFzEvent digitize_calo_muon_rich(
    std::int32_t run,
    std::int32_t event,
    std::span<const SimShowerParticle> particles,
    std::uint64_t seed = 42);

}  // namespace delphipp::sim
