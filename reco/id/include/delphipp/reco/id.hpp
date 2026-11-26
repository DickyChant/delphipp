// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Inner Detector (ID) reconstruction
// (simana_v94d/src/car/idana.car: 24-sector jet chamber drift hits + 5 trigger
// MWPC layers -> local (r, r*phi) track segments).

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::id {

inline constexpr std::int32_t kDetectorId = 4;  // ID detector code
inline constexpr std::int32_t kNumSectors = 24;

struct IdJetHit {
  std::int32_t sector = 1;  // 1..24
  std::int32_t wire = 1;    // 1..24
  float r_cm = 15.0f;       // 11.8 .. 22.3 cm
  float drift_time_ns = 0.0f;
  int lr_sign = +1;         // Left/right ambiguity sign (+1 or -1)
};

struct IdSegment {
  std::int32_t sector = 1;
  float r_ref_cm = 17.0f;
  float phi_rad = 0.0f;
  float dphi_dr = 0.0f;
  std::int32_t num_hits = 0;
};

[[nodiscard]] std::vector<IdJetHit> unpack_raw_blocklet(
    std::span<const std::int32_t> words);

[[nodiscard]] std::vector<IdSegment> reconstruct_segments(
    std::span<const IdJetHit> hits,
    float v_drift_cm_per_ns = 0.0040f);

std::vector<TanagraId> write_ter_banks(
    TanagraStore& store,
    std::span<const IdSegment> segments);

}  // namespace delphipp::reco::id
