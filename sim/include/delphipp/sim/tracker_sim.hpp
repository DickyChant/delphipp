// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELSIM tracker digitization & MC truth bank builder
// (simana_v94d/src/car/vdsim.car, idsim.car, tpcsim.car, odsim.car,
// fcasim.car, fcbsim.car, stwana.car, delsim36.car):
//   - Produces relocated TOP -> RAW -> RVD/RID/RTPC/ROD/RFCA/RFCB -> RAWD
//     (plus STWR) ZEBRA banks decoding 1-to-1 through all 6 tracker unpackers.
//   - Produces DELSIM SIM_ -> SP -> ST truth records.

#pragma once

#include <cstdint>

#include "delphipp/io/fz_zebra_parser.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/sim/sim_core.hpp"

namespace delphipp::sim {

// Digitize all particles in `ev` across VD, ID, STW, TPC, OD, FCA, and FCB,
// returning a relocated ParsedFzEvent with TOP -> RAW -> PMB -> RAWD (+ STWR)
// banks ready for EventReconstructionPipeline and raw_blocklet unpackers.
[[nodiscard]] io::ParsedFzEvent digitize_trackers(
    const SimEvent& ev,
    const MagneticField& bfield,
    std::uint64_t seed = 1ULL);

// Build a DELSIM MC truth event containing relocated SIM_ -> SP -> ST banks.
[[nodiscard]] io::ParsedFzEvent build_sim_truth_record(const SimEvent& ev);

}  // namespace delphipp::sim
