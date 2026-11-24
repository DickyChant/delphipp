// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI Forward Muon Chambers (MUF) reconstruction
// (simana_v94d/src/car/mufana.car: EQTSLW, EQDLDF, EQWPOS, EQDRIF, EQSPAF).

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace delphipp::reco::muf {

inline constexpr std::int32_t kDetectorId = 30;  // MUF detector code
inline constexpr std::size_t kNumQuadrants = 16;

// Mirrors COMMON /EQCONS/ initialized in EQINI (mufana.car lines 89-118)
struct MufConfig {
  float tsum_window_ns = 50.0f;                                   // TSUMTE
  std::array<float, 2> ctslw0 = {-7.680f, -5.12f};                // long, short
  std::array<float, 2> ctslw1 = {0.964f, 0.59f};
  std::array<float, 2> ctslw2 = {2.643f, 3.17f};
  std::array<float, 2> ctslw3 = {-0.367f, -0.33f};
  float sigma_drift_cm = 0.5f;                                    // DDRIFT (RMSDRI)
  float sigma_delay_line_cm = 0.5f;                               // DXONL2 (RMSDLI)
  float cut_lr = 3.0f;                                            // CUTLR
  float bmn_ltd = 24.9f;                                          // BMNLTD
  float prepulse_cut_ns = -25.0f;                                 // PREPUL
  float flight_time_ns = 20.0f;                                   // TFLI in EQSPAF
  bool is_simulation = false;                                     // ISIMUL != 0
};

// Subset of /EQGDDB/ + /EQCDDB/ needed for local hit reconstruction
struct MufCalibration {
  float vanode_cm_per_ns = 25.0f;  // VANODE
  float wire_pitch_cm = 18.0f;     // STPWBW
  // [quad 0..15] -> QGEOM(1..6, quad): [2]=x0, [3]=y0, [4]=angle_x, [5]=angle_y
  std::array<std::array<float, 6>, kNumQuadrants> qgeom{};
  // [quad 0..15][layer 0..1] -> SVDRIF(layer, quad) in cm/ns
  std::array<std::array<float, 2>, kNumQuadrants> svdrif{};
  // [quad 0..15][layer 0..1] -> SX0PAR(layer, quad) in cm
  std::array<std::array<float, 2>, kNumQuadrants> sx0par{};
};

struct MufRawTriplet {
  std::int32_t quadrant = 1;  // 1..16 (NQUAD)
  std::int32_t chamber = 1;   // 1..22 (NCHINL)
  std::int32_t layer = 1;     // 1=X, 2=Y (IHV)
  float ta_ns = 0.0f;         // anode time (RTA)
  float tl1_ns = -1.0f;       // delay-line end 1 time (<0 if missing)
  float tl2_ns = -1.0f;       // delay-line end 2 time (<0 if missing)
  float tca_ns = 0.0f;        // anode cable delay
  float tc1_ns = 0.0f;        // end-1 cable delay
  float tc2_ns = 0.0f;        // end-2 cable delay
  float vline_cm_per_ns = 0.0f;
  float chamber_len_cm = 0.0f;
  float time_sum_ns = 0.0f;   // TEST
};

struct MufSpaceHit {
  std::int32_t quadrant = 0;
  std::int32_t chamber = 0;
  std::int32_t layer = 0;
  std::int32_t encoded_id = 0;  // quad*1000 + chamber*10 + layer
  float drift_dist_cm = 0.0f;
  float x_on_line_cm = 0.0f;
  float x_wire_cm = 0.0f;
  float sigma_drift_cm = 0.0f;
  float sigma_line_cm = 0.0f;
};

struct DelayLineResult {
  float x_on_line_cm = -100.0f;
  float sigma_cm = 0.0f;
};

// Port of SUBROUTINE EQTSLW
[[nodiscard]] float time_slew(int end_1or2,
                              float chamber_len_cm,
                              float vline_cm_per_ns,
                              float time_ns,
                              const MufConfig& cfg = {}) noexcept;

// Port of SUBROUTINE EQDLDF
[[nodiscard]] DelayLineResult delay_line_distance(const MufRawTriplet& hit,
                                                  float vanode_cm_per_ns,
                                                  const MufConfig& cfg = {}) noexcept;

// Port of SUBROUTINE EQWPOS
[[nodiscard]] float wire_position(std::int32_t quadrant_1based,
                                  std::int32_t chamber_1based,
                                  std::int32_t layer_1or2,
                                  float x_on_line_cm,
                                  const MufCalibration& calib) noexcept;

// Port of SUBROUTINE EQDRIF
[[nodiscard]] float drift_distance(std::int32_t quadrant_1based,
                                   std::int32_t layer_1or2,
                                   float drift_time_ns,
                                   const MufCalibration& calib) noexcept;

// Port of SUBROUTINE EQSPAF
[[nodiscard]] std::vector<MufSpaceHit> reconstruct_hits(
    std::span<const MufRawTriplet> triplets,
    const MufCalibration& calib,
    const MufConfig& cfg = {});

}  // namespace delphipp::reco::muf
