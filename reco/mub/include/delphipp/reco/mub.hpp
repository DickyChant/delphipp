// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Barrel Muon Chambers (MUB) reconstruction
// (simana_v94d/src/car/mubana.car: EMCHAN, EMCHID, EMSLPB, EMLTIM, EMDFTT).

#pragma once

#include <array>
#include <cstdint>
#include <optional>

namespace delphipp::reco::mub {

inline constexpr std::int32_t kDetectorId = 14;    // MUB detector code
inline constexpr std::int32_t kNumChambers = 1372; // 48*28 + 2*14

struct ChamberLocation {
  std::int32_t sector = 0;  // 1..24 End A, 25..48 End C, 49..50 Leg
  std::int32_t layer = 0;   // 1..3 Inner, 4..5 Outer, 6..7 Peripheral
  std::int32_t pos = 0;     // 1..(3,4,5)
};

// Calibration constants for one LTD module (CONLTD(1..9, iltd))
struct LtdCalibration {
  float coarse_step_ns = 0.0f;            // CONLTD(6)
  std::array<float, 3> fine_step_ns{};    // CONLTD(7..9)
};

// Calibration constants for one MUB chamber (from /EMEDDB/ & /EMCOMB/)
struct MubChamberCalib {
  float length_cm = 365.0f;                  // CHBLEN
  float half_width_cm = 10.4f;               // CWID2
  std::array<float, 3> tzero_ns{};           // CTZERO(ich, 1..3): anode, near, far
  std::array<float, 5> delay_poly{};         // CDELAY(ich, 1, 1..5)
  std::array<float, 4> anode_poly{};         // CANODE(ich, 1..4)
};

struct MubConfig {
  std::array<float, 7> tof_ns{};   // EMTOFL(1..7)
  float vanode_cm_per_ns = 25.0f;  // VANODE
  float veldrf_cm_per_ns = 0.005f; // VELDRF
  float atmosp_mbar = 980.0f;      // ATMOSP (/EXMACH/)
  float tmax_anode_ns = 4000.0f;   // TMAX
  float tmax_dl_ns = 6000.0f;      // TMAX1
  float cosmic_bcross_ns = 0.0f;   // BCROSS
  bool cosmic_mode = false;        // COSFLG
  bool is_simulation = false;      // IFSIM == 1
};

struct MubTripletResult {
  float z_local_cm = 0.0f;   // Cargo convention: 0 at far end, +ve to near end
  float drift_cm = 0.0f;     // DDRFT
  float anode_time_ns = 0.0f;// TDD
  float time_sum_ns = 0.0f;  // TSUM
};

// Port of SUBROUTINE EMCHAN: (sector, layer, pos) -> 1..1372 or nullopt if invalid
[[nodiscard]] std::optional<std::int32_t> chamber_index(
    std::int32_t sector, std::int32_t layer, std::int32_t pos) noexcept;

// Port of SUBROUTINE EMSLPB: inverse of EMCHAN (1..1372 -> ChamberLocation)
[[nodiscard]] std::optional<ChamberLocation> chamber_location(
    std::int32_t index_1based) noexcept;

// Port of SUBROUTINE EMCHID: (imub 1..3, isect 1..24, ihemb 1..3, icham) -> 1..1372
[[nodiscard]] std::int32_t chamber_id(std::int32_t imub,
                                      std::int32_t isect,
                                      std::int32_t ihemb,
                                      std::int32_t icham) noexcept;

// Port of SUBROUTINE EMLTIM: raw LTD word -> calibrated time in ns
[[nodiscard]] float ltd_time_ns(std::int32_t raw_ltd_word,
                                const LtdCalibration& calib) noexcept;

// Port of SUBROUTINE EMDFTT: triplet (anode, near, far) -> local z & drift distance
[[nodiscard]] std::optional<MubTripletResult> reconstruct_triplet(
    std::int32_t layer_1based,
    std::array<float, 3> times_ns,
    const MubChamberCalib& chamber,
    const MubConfig& cfg = {}) noexcept;

}  // namespace delphipp::reco::mub
