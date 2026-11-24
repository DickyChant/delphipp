// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// C++20 port of DELPHI Forward Hodoscope (HOF) reconstruction (simana_v94d/src/car/hofana.car).
// Preserves exact bit-unpacking, coordinate, covariance, timing, and bunch-ID semantics.

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/zebra_store.hpp"

namespace delphipp::reco::hof {

inline constexpr std::int32_t kDetectorId = 29;       // HOF detector code
inline constexpr std::int32_t kModuleRaw = 29;        // Real-data raw module ID
inline constexpr std::int32_t kModuleSim = 22;        // Simulated raw module ID
inline constexpr std::int32_t kModuleTiming = 30;     // EQRAW timing module ID
inline constexpr std::size_t kGeoDdbWords = 1120;     // 2 endcaps * 4 quads * 28 counters * 5 floats
inline constexpr std::size_t kCalDdbWords = 278;      // 256 pad channels + 2 * 11 timing channels
inline constexpr std::size_t kMaxHits = 250;

// Mirrors COMMON /EGSTAT/ NEVEG, NSIGN, NRAWBR, NSGQUA(8), NSGFUL(12)
struct HofStats {
  std::int32_t nevegs = 0;
  std::int32_t nsign = 0;
  std::int32_t nrawbr = 0;
  std::array<std::int32_t, 8> nsgqua{};
  std::array<std::int32_t, 12> nsgful{};

  void reset_run() noexcept {
    nevegs = 0;
    nsign = 0;
    nrawbr = 0;
    nsgqua.fill(0);
    nsgful.fill(0);
  }
  void begin_event() noexcept {
    nsign = 0;
    ++nevegs;
  }
};

// Mirrors COMMON /EGGDDB/ GEODDB(1120) and /EGCDDB/ CALDDB(278)
struct HofCalibration {
  std::array<float, kGeoDdbWords> geoddb{};
  std::array<std::int32_t, kCalDdbWords> ichan{};
};

// One 8-word TD hit entry (DATTD(1..8, ihit))
struct HofTdHit {
  std::array<float, 8> words{};

  [[nodiscard]] bool is_counter_hit() const noexcept { return words[0] != 0.0f; }
  [[nodiscard]] float x() const noexcept { return words[1]; }
  [[nodiscard]] float y() const noexcept { return words[2]; }
  [[nodiscard]] float z() const noexcept { return words[3]; }
  [[nodiscard]] float dx() const noexcept { return words[4]; }
  [[nodiscard]] float dy() const noexcept { return words[5]; }
  [[nodiscard]] float dz() const noexcept { return words[6]; }
  [[nodiscard]] float id_word() const noexcept { return words[7]; }
};

// One 24-word TE bank payload (DATTE(1..24) / IDATTE(1..24))
struct HofTeBank {
  std::array<std::int32_t, 24> raw_words{};

  void set_int(std::size_t one_based_idx, std::int32_t val) noexcept;
  void set_float(std::size_t one_based_idx, float val) noexcept;
  [[nodiscard]] std::int32_t get_int(std::size_t one_based_idx) const noexcept;
  [[nodiscard]] float get_float(std::size_t one_based_idx) const noexcept;
};

struct HofTdResult {
  std::array<std::int32_t, 9> header_raw{};
  std::vector<HofTdHit> hits;
  std::int32_t counter_hits = 0;  // NHITS (counter hits only; hits.size() == NHIT)
};

// Bit-extraction helpers matching Fortran IBITS/BTEST (0-based) and CERNLIB JBIT/JBYT (1-based)
[[nodiscard]] constexpr std::uint32_t ibits(std::int32_t word, int pos, int len) noexcept {
  const auto u = static_cast<std::uint32_t>(word);
  const std::uint32_t mask = (len >= 32) ? 0xFFFFFFFFu : ((1u << len) - 1u);
  return (u >> pos) & mask;
}
[[nodiscard]] constexpr bool btest(std::int32_t word, int pos) noexcept {
  return ((static_cast<std::uint32_t>(word) >> pos) & 1u) != 0u;
}
[[nodiscard]] constexpr std::uint32_t jbit(std::int32_t word, int loc_1based) noexcept {
  return ibits(word, loc_1based - 1, 1);
}
[[nodiscard]] constexpr std::uint32_t jbyt(std::int32_t word, int loc_1based, int len) noexcept {
  return ibits(word, loc_1based - 1, len);
}

// Port of SUBROUTINE EGTD(NHITS): unpack EGRAW (lraw) and EQRAW (mraw) into TD hits.
[[nodiscard]] HofTdResult unpack_td(const ZebraStoreView& store,
                                    std::int32_t lraw,
                                    std::int32_t mraw,
                                    const HofCalibration& calib,
                                    HofStats& stats);

// Port of SUBROUTINE EGTE(NHITS): build TE banks from counter hits and EQRAW timing.
[[nodiscard]] std::vector<HofTeBank> build_te(const ZebraStoreView& store,
                                              std::int32_t eqraw_link,
                                              const HofTdResult& td,
                                              const HofCalibration& calib);

// Port of SUBROUTINE EGNBUN(NBUNCH): determine LEP bunch number from HOF TE times.
[[nodiscard]] std::int32_t identify_bunch(std::span<const HofTeBank> tes,
                                          std::int32_t nbulep,
                                          float bunint_ns);

}  // namespace delphipp::reco::hof
