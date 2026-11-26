// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Unified pure-C++20 DELPHI event reconstruction pipeline orchestrating
// ZEBRA store population, subdetector raw blocklet decoding (18 detectors),
// TANAGRA TDR/TER/TKR banking, ambiguity resolution, Kalman refit, vertexing,
// and Short-DST construction.

#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>

#include "delphipp/dst/short_dst.hpp"
#include "delphipp/io/fz_zebra_parser.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/tanagra_store.hpp"

namespace delphipp::pipeline {

// Map DELPHI numeric detector ID (e.g. 5 -> "TPC", 29 -> "HOF") to short name.
[[nodiscard]] std::string detector_name(std::int32_t detector_id);

struct EventRecoOutput {
  std::int32_t run = 0;
  std::int32_t event = 0;
  std::string top_bank;
  std::size_t zebra_banks = 0;
  std::size_t raw_blocklets = 0;
  std::int32_t hof_hits = 0;
  std::map<std::string, std::size_t> detector_blocklet_counts;
  TanagraStore tanagra;
  dst::DstEvent dst_event;
  std::size_t tracks_before_ambiguity = 0;
  std::size_t tracks_after_ambiguity = 0;
  std::size_t kalman_refits = 0;
};

class EventReconstructionPipeline {
 public:
  explicit EventReconstructionPipeline(const MagneticField& bfield,
                                       float abst0_val = 0.0f) noexcept;

  [[nodiscard]] EventRecoOutput reconstruct_event(
      const io::ParsedFzEvent& ev,
      std::int32_t fallback_run = 47014,
      std::int32_t fallback_evt = 1) const;

 private:
  const MagneticField& bfield_;
  float abst0_val_ = 0.0f;
};

}  // namespace delphipp::pipeline
