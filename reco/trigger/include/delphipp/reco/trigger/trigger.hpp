// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 port of DELPHI T1-T4 trigger simulation and decoding
// (simana_v94d/src/car/trigger.car & fcabtrig.car):
//  - T1 hardware subtrigger coincidences (ID, OD, TOF, FCA, FCB, HPC, EMF, MUB, SAT)
//  - FCA-FCB forward/backward 3-projection (X, U, V) T1/T2 trigger (VA33/VA40)
//  - T2 TPC / HPC / HCAL / FCB confirmations
//  - T3 / T4 physics stream classification (HadronicZ, Dilepton, Bhabha, GammaGamma)

#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace delphipp::reco::trigger {

// FCA-FCB T1/T2 forward track trigger inputs (fcabtrig.car)
struct FcabEndcapInput {
  // Number of hit planes (0..4) in each FCB coordinate projection (V, U, X)
  std::array<std::int32_t, 3> fcb_plane_hits{};
  // Group numbers (0..63) after T1B per projection (V, U, X) for T2 sum check
  std::array<std::vector<std::int32_t>, 3> t1b_groups{};
  std::int32_t fca_hit_groups = 0;
};

struct FcabTriggerInput {
  FcabEndcapInput forward{};   // +z endcap
  FcabEndcapInput backward{};  // -z endcap
};

struct FcabTriggerResult {
  std::int32_t t1_forward = 0;       // IFFW (0..3)
  std::int32_t t1_backward = 0;      // IFBW (0..3)
  std::int32_t t2_forward = 0;       // ISFW (0..3), count of 93 <= sum <= 100
  std::int32_t t2_backward = 0;      // ISBW (0..3)
  std::int32_t n_3coord_forward = 0; // Projections meeting majority on +z
  std::int32_t n_3coord_backward = 0;
};

// Subdetector trigger primitives feeding DELPHI T1, T2, T3, T4
struct TriggerDetectorInput {
  // T1 tracking & timing chambers
  std::int32_t id_hits = 0;
  std::int32_t od_hits = 0;
  std::int32_t tof_hits = 0;
  bool tof_back_to_back = false;

  // Forward chambers FCA/FCB (if populated, overrides fca/fcb_coinc_* below)
  FcabTriggerInput fcab{};
  bool fca_fcb_fw = false;
  bool fca_fcb_bw = false;
  bool fcab_t2_fw = false;
  bool fcab_t2_bw = false;

  // Barrel & forward calorimetry + HOF + luminosity
  std::int32_t hpc_clusters = 0;
  std::int32_t hpc_high_clusters = 0;
  float hpc_energy_gev = 0.0f;
  bool hpc_back_to_back = false;

  bool emf_low_fw = false;
  bool emf_low_bw = false;
  bool emf_high_fw = false;
  bool emf_high_bw = false;
  float emf_energy_fw_gev = 0.0f;
  float emf_energy_bw_gev = 0.0f;

  std::int32_t hof_hits_fw = 0;
  std::int32_t hof_hits_bw = 0;
  bool sat_bhabha_coincidence = false;

  // Muon chambers
  std::int32_t mub_hits = 0;
  bool mub_back_to_back = false;
  std::int32_t muf_hits = 0;

  // T2 confirmation detectors (TPC & HCAL)
  std::int32_t tpc_tracks = 0;
  bool tpc_back_to_back = false;
  std::int32_t hcal_clusters = 0;
  float hcal_energy_gev = 0.0f;

  // T3 / T4 online/offline reconstructed summary
  std::int32_t charged_multiplicity = 0;
  float visible_energy_gev = 0.0f;
  float max_track_p_gev = 0.0f;
};

// Optional T1 subtrigger mask (true = enabled)
struct T1Mask {
  bool id_od = true;
  bool tof_od = true;
  bool tof_bb = true;
  bool hpc_od = true;
  bool hpc_tof = true;
  bool hpc_bb = true;
  bool hpc_maj = true;
  bool fwd_maj = true;
  bool emf_bb = true;
  bool mub_coinc = true;
  bool sat_bhabha = true;
};

struct T1Result {
  bool id_od = false;       // ID >= 1 * OD >= 2 (IDOD)
  bool tof_od = false;      // TOF >= 1 * OD >= 1 (SCOD)
  bool tof_bb = false;      // TOF >= 2 back-to-back (TOBB)
  bool hpc_od = false;      // HPC >= 1 * OD >= 1 (HPOD)
  bool hpc_tof = false;     // HPC >= 1 * TOF >= 1 (THMJ)
  bool hpc_bb = false;      // HPC >= 2 back-to-back (HPBB)
  bool hpc_maj = false;     // HPC majority (HPM1 / HPM2)
  bool fwd_maj = false;     // Forward majority (FWMJ: FCA*FCB*EMF/HOF)
  bool emf_bb = false;      // EMF forward-backward coincidence (FEBH / FEHI)
  bool mub_coinc = false;   // MUB * (ID || OD || TOF) or MUB back-to-back
  bool sat_bhabha = false;  // SAT Bhabha arm coincidence (SABH / SASV)
  bool t1_accept = false;
};

struct T2Result {
  bool tpc_confirm = false;   // TPC T2 tracks (TPLO / TPBB / TPGL)
  bool hpc_confirm = false;   // HPC T2 energy/cluster confirmation
  bool hcal_confirm = false;  // HCAL T2 energy/cluster confirmation
  bool fcab_confirm = false;  // FCA-FCB T2S collinearity sum confirmation
  bool mub_confirm = false;   // MUB T2 track coincidence
  bool emf_sat_confirm = false;
  bool t2_accept = false;
};

enum class PhysicsStream : std::uint32_t {
  None       = 0,
  HadronicZ  = 1u << 0,
  Dilepton   = 1u << 1,
  Bhabha     = 1u << 2,
  GammaGamma = 1u << 3,
};

[[nodiscard]] constexpr PhysicsStream operator|(PhysicsStream a,
                                                PhysicsStream b) noexcept {
  return static_cast<PhysicsStream>(static_cast<std::uint32_t>(a) |
                                    static_cast<std::uint32_t>(b));
}

[[nodiscard]] constexpr PhysicsStream operator&(PhysicsStream a,
                                                PhysicsStream b) noexcept {
  return static_cast<PhysicsStream>(static_cast<std::uint32_t>(a) &
                                    static_cast<std::uint32_t>(b));
}

[[nodiscard]] constexpr bool has_stream(PhysicsStream mask,
                                        PhysicsStream stream) noexcept {
  return (static_cast<std::uint32_t>(mask) &
          static_cast<std::uint32_t>(stream)) != 0u;
}

struct TriggerDecision {
  FcabTriggerResult fcab{};
  T1Result t1{};
  T2Result t2{};
  bool t3_accept = false;
  bool t4_accept = false;
  PhysicsStream streams = PhysicsStream::None;
  PhysicsStream primary_stream = PhysicsStream::None;
};

// Evaluate FCA-FCB T1/T2 trigger (fcabtrig.car VATR1 / VAEND).
[[nodiscard]] FcabTriggerResult evaluate_fcab_trigger(
    const FcabTriggerInput& input) noexcept;

// Evaluate full DELPHI T1 -> T2 -> T3 -> T4 trigger chain.
[[nodiscard]] TriggerDecision evaluate_trigger(
    const TriggerDetectorInput& input,
    const T1Mask& mask = {}) noexcept;

}  // namespace delphipp::reco::trigger
