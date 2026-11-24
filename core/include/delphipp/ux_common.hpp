// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Binary-compatible C++20 mappings for core DELPHI UX / ZEBRA COMMON blocks,
// derived from simana_v94d/src/car/ux26.car (lines 432-508).

#pragma once

#include <cstddef>
#include <cstdint>

namespace delphipp::fortran {

inline constexpr std::size_t kDefaultNSizeQ = 2'000'000;
inline constexpr std::size_t kDefaultNSizWs = 500'000;
inline constexpr std::size_t kDefaultNWTans = 1'000'000;
inline constexpr std::size_t kNumDet = 20;

// COMMON /QUEST/ IQUEST(100)
struct QuestCommon {
  std::int32_t iquest[100];
};
static_assert(sizeof(QuestCommon) == 100 * sizeof(std::int32_t));

// COMMON /UXCDIV/ IUXSTO, IUXUSE, IUXRAW, IUXSIM, IUXINT,
//                 IUXDST, IUXRUN, IUXTIT, IUXCAL, IUXDET, IUXWIP
struct UxDivCommon {
  std::int32_t iuxsto;
  std::int32_t iuxuse;
  std::int32_t iuxraw;
  std::int32_t iuxsim;
  std::int32_t iuxint;
  std::int32_t iuxdst;
  std::int32_t iuxrun;
  std::int32_t iuxtit;
  std::int32_t iuxcal;
  std::int32_t iuxdet;
  std::int32_t iuxwip;
};
static_assert(sizeof(UxDivCommon) == 11 * sizeof(std::int32_t));

// COMMON /UXLINK/ LTEMP(2), LRTOP, LSTOP, LTTOP, LITOP, LRTEMP,
//                 LRWTMP, LRAWUX, LBKTOP, LORTOP, LRTINT, LDTOP
struct UxLinkCommon {
  std::int32_t ltemp[2];
  std::int32_t lrtop;
  std::int32_t lstop;
  std::int32_t lttop;
  std::int32_t litop;
  std::int32_t lrtemp;
  std::int32_t lrwtmp;
  std::int32_t lrawux;
  std::int32_t lbktop;
  std::int32_t lortop;
  std::int32_t lrtint;
  std::int32_t ldtop;
};
static_assert(sizeof(UxLinkCommon) == 13 * sizeof(std::int32_t));

// COMMON /UXCEVT/ NEXP, NROLD, IEVOLD, NRUN, IEVT, NFILSQ,
//                 IEVDAT, IEVTIM, NEVT, NEVTOT, NCPT, NFIOLD,
//                 NFILL, NVERDS, IPRDAT, IPRTIM
struct UxEvtCommon {
  std::int32_t nexp;
  std::int32_t nrold;
  std::int32_t ievold;
  std::int32_t nrun;
  std::int32_t ievt;
  std::int32_t nfilsq;
  std::int32_t ievdat;
  std::int32_t ievtim;
  std::int32_t nevt;
  std::int32_t nevtot;
  std::int32_t ncpt;
  std::int32_t nfiold;
  std::int32_t nfill;
  std::int32_t nverds;
  std::int32_t iprdat;
  std::int32_t iprtim;
};
static_assert(sizeof(UxEvtCommon) == 16 * sizeof(std::int32_t));

// COMMON /UXCOM/ UXFENC(10), LUXSTR(30), LUXREF(30)
// With EQUIVALENCE (LQ(1), LUXSTR(1)), (Q(1), IQ(1), LQ(9)).
// Note: In Fortran, LQ(1) starts at word index 10 of /UXCOM/ (immediately after
// UXFENC(10)), and IQ(1)/Q(1) start at LQ(9), i.e., word index 18 of /UXCOM/.
inline constexpr std::size_t kUxFenceWords = 10;
inline constexpr std::size_t kLqToIqOffset = 8;  // LQ(9) == IQ(1) == Q(1)

}  // namespace delphipp::fortran
