// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/hof.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <numbers>

namespace delphipp::reco::hof {

namespace {

// Bit masks from hofana.car
constexpr std::int32_t kIgdscr = 8 + 16 + 32 + 64 + 256;          // EGTD
constexpr std::int32_t kIpdscr = 2 + 4 + 8;                       // EGTD
constexpr std::int32_t kIdescr = 4 + 8 + 16 + 32 + 64 + 16384;    // EGTE
constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;

constexpr std::int32_t kLtdTab[5] = {1, 2, 3, 4, 5};

void set_header_float(std::array<std::int32_t, 9>& hdr, std::size_t one_based, float val) noexcept {
  hdr[one_based - 1] = std::bit_cast<std::int32_t>(val);
}

}  // namespace

void HofTeBank::set_int(std::size_t one_based_idx, std::int32_t val) noexcept {
  raw_words[one_based_idx - 1] = val;
}

void HofTeBank::set_float(std::size_t one_based_idx, float val) noexcept {
  raw_words[one_based_idx - 1] = std::bit_cast<std::int32_t>(val);
}

std::int32_t HofTeBank::get_int(std::size_t one_based_idx) const noexcept {
  return raw_words[one_based_idx - 1];
}

float HofTeBank::get_float(std::size_t one_based_idx) const noexcept {
  return std::bit_cast<float>(raw_words[one_based_idx - 1]);
}

HofTdResult unpack_td(const ZebraStoreView& store,
                      std::int32_t lraw,
                      std::int32_t mraw,
                      const HofCalibration& calib,
                      HofStats& stats) {
  HofTdResult out;
  out.header_raw[0] = kDetectorId;
  out.header_raw[1] = 0;
  out.header_raw[2] = kIgdscr;
  out.header_raw[4] = 8;
  out.header_raw[5] = 0;
  set_header_float(out.header_raw, 7, 0.0f);
  set_header_float(out.header_raw, 8, 0.0f);
  set_header_float(out.header_raw, 9, 0.5f);

  // 1. Unpack counter hits from EGRAW linear bank chain
  std::int32_t loc = lraw;
  while (loc != 0) {
    if (store.iq(loc - 1) <= 5) {
      break;
    }
    const auto module = static_cast<std::int32_t>(ibits(store.iq(loc + 2), 8, 8));
    const auto numcap = static_cast<std::int32_t>(ibits(store.iq(loc + 2), 16, 4));

    if ((module == kModuleRaw || module == kModuleSim) && (numcap == 1 || numcap == 2)) {
      for (std::int32_t iquadr = 1; iquadr <= 4; ++iquadr) {
        const std::int32_t infhof = store.iq(loc + 4 + iquadr);
        if (infhof == 0) {
          continue;
        }
        const std::size_t ipnter =
            static_cast<std::size_t>((numcap - 1) * 128 + (iquadr - 1) * 32);
        for (int ic = 0; ic < 32; ++ic) {
          if (!btest(infhof, ic)) {
            continue;
          }
          std::int32_t icount = 0;
          if (module == kModuleSim) {
            icount = ic + 1;
          } else if (module == kModuleRaw) {
            icount = calib.ichan[ipnter + static_cast<std::size_t>(ic)];
          }
          if (icount < 1 || icount > 28) {
            continue;
          }
          if (out.hits.size() >= kMaxHits) {
            break;
          }
          HofTdHit hit;
          hit.words[0] = static_cast<float>(kIpdscr);
          const std::size_t ipnt =
              static_cast<std::size_t>(((numcap - 1) * 112 + (iquadr - 1) * 28 + icount - 1) * 5);
          for (std::size_t k = 0; k < 5; ++k) {
            hit.words[1 + k] = calib.geoddb[ipnt + k];
          }
          hit.words[6] = 0.5f;
          hit.words[7] = static_cast<float>(numcap * 1000 + iquadr * 100 + icount);
          out.hits.push_back(hit);

          const std::size_t iwrk = static_cast<std::size_t>((numcap - 1) * 4 + (iquadr - 1));
          stats.nsgqua[iwrk] += 1;
        }
      }
    }
    loc = store.lq(loc);
  }

  out.counter_hits = static_cast<std::int32_t>(out.hits.size());

  // 2. Unpack timing hits from EQRAW (mraw) linear bank chain
  if (mraw > 0 && store.iq(mraw - 1) > 5) {
    // Note: hofana.car line 549 reads NEWPNT from LRAW+2 if LRAW > 0
    const auto newpnt = (lraw > 0)
                            ? static_cast<std::int32_t>(jbit(store.iq(lraw + 2), 1))
                            : static_cast<std::int32_t>(jbit(store.iq(mraw + 2), 1));
    loc = mraw;
    while (loc != 0) {
      if (store.iq(loc - 1) <= 5) {
        break;
      }
      const auto nendcp = static_cast<std::int32_t>(jbyt(store.iq(loc + 2), 17, 4));
      const bool valid_hdr = (nendcp == 1 || nendcp == 2) &&
                             (jbyt(store.iq(loc + 2), 21, 12) == 0) &&
                             (jbyt(store.iq(loc + 2), 9, 8) == kModuleTiming);
      if (valid_hdr) {
        for (int ll = 0; ll < 5; ++ll) {
          const std::int32_t iii = kLtdTab[ll];
          const std::int32_t ites_hdr = store.iq(loc + iii + 4);
          if (ites_hdr == 0) {
            continue;
          }
          const auto iwd1 = static_cast<std::int32_t>(jbyt(ites_hdr, 1, 17)) + newpnt;
          const auto nh = static_cast<std::int32_t>(jbyt(ites_hdr, 18, 10));
          const auto idrop = static_cast<std::int32_t>(jbyt(ites_hdr, 28, 5));
          if (idrop != 0 || nh == 0) {
            continue;
          }
          const auto islot = static_cast<std::int32_t>(jbyt(store.iq(loc + iwd1 + 3), 17, 5));
          const std::int32_t locw = loc + iwd1 + 5;
          const std::size_t ipoint = static_cast<std::size_t>(256 + 11 * (nendcp - 1));

          for (std::int32_t kl = 1; kl <= nh; ++kl) {
            const std::int32_t ites = store.iq(locw + kl);
            const auto nchan = static_cast<std::int32_t>(jbyt(ites, 15, 6));
            for (int ipp = 1; ipp <= 4; ++ipp) {
              const std::int32_t isldb = calib.ichan[ipoint + static_cast<std::size_t>(ipp)];
              const std::int32_t ichdb = calib.ichan[ipoint + static_cast<std::size_t>(ipp + 5)];
              if (islot == isldb && nchan == ichdb && out.hits.size() < kMaxHits) {
                const auto ihoft = static_cast<std::int32_t>(jbyt(ites, 1, 14));
                HofTdHit thit;
                thit.words[7] = static_cast<float>(ihoft * 100 + nendcp * 10 + ipp);
                out.hits.push_back(thit);
              }
            }
          }
        }
      }
      loc = store.lq(loc);
    }
  }

  out.header_raw[3] = static_cast<std::int32_t>(out.hits.size());

  const std::int32_t nhitl = std::min<std::int32_t>(out.counter_hits, 11);
  stats.nsgful[static_cast<std::size_t>(nhitl)] += 1;

  return out;
}

std::vector<HofTeBank> build_te(const ZebraStoreView& store,
                                std::int32_t eqraw_link,
                                const HofTdResult& td,
                                const HofCalibration& calib) {
  std::vector<HofTeBank> tes;
  const auto n_counter = static_cast<std::size_t>(
      std::min<std::int32_t>(td.counter_hits, static_cast<std::int32_t>(td.hits.size())));
  tes.reserve(n_counter);

  for (std::size_t i = 0; i < n_counter; ++i) {
    const HofTdHit& hit = td.hits[i];
    HofTeBank te;

    const auto iduser = static_cast<std::int32_t>(hit.id_word());
    const auto modsub = static_cast<std::int32_t>(hit.id_word() / 100.0f);

    te.set_int(1, kDetectorId);
    te.set_int(2, modsub);
    te.set_int(3, 0);
    te.set_int(4, kIdescr);
    te.set_int(5, 21);
    te.set_int(6, 3);
    te.set_int(7, 0);
    te.set_float(8, 0.0f);
    te.set_float(9, 1.0f);

    // Coordinates (x, y, z) and spherical angles (theta, phi)
    const float x = hit.x();
    const float y = hit.y();
    const float z = hit.z();
    te.set_float(10, x);
    te.set_float(11, y);
    te.set_float(12, z);
    const float r = std::hypot(x, y);
    te.set_float(13, std::atan2(r, z));
    float phi = std::atan2(y, x);
    if (phi < 0.0f) {
      phi += kTwoPi;
    }
    te.set_float(14, phi);
    te.set_float(15, 0.0f);

    // Timing lookup in EQRAW (take first matching timing hit for this endcap+quadrant)
    float timing = 0.0f;
    const std::int32_t numcap = iduser / 1000;
    const std::int32_t iquadr = (iduser % 1000) / 100;
    if ((numcap == 1 || numcap == 2) && (iquadr >= 1 && iquadr <= 4) &&
        eqraw_link > 0 && store.iq(eqraw_link - 1) > 5) {
      const auto newpnt = static_cast<std::int32_t>(jbit(store.iq(eqraw_link + 2), 1));
      std::int32_t loc = eqraw_link;
      bool found = false;
      while (loc != 0 && !found) {
        if (store.iq(loc - 1) <= 5) {
          break;
        }
        const auto nendcp = static_cast<std::int32_t>(jbyt(store.iq(loc + 2), 17, 4));
        const bool match_hdr = (nendcp == numcap) &&
                               (jbyt(store.iq(loc + 2), 21, 12) == 0) &&
                               (jbyt(store.iq(loc + 2), 9, 8) == kModuleTiming);
        if (match_hdr) {
          const std::size_t ipoint = static_cast<std::size_t>(256 + 11 * (numcap - 1));
          const std::int32_t isldb = calib.ichan[ipoint + static_cast<std::size_t>(iquadr)];
          const std::int32_t ichdb = calib.ichan[ipoint + static_cast<std::size_t>(iquadr + 5)];
          for (int ll = 0; ll < 5 && !found; ++ll) {
            const std::int32_t iii = kLtdTab[ll];
            const std::int32_t ites_hdr = store.iq(loc + iii + 4);
            if (ites_hdr == 0) {
              continue;
            }
            const auto iwd1 = static_cast<std::int32_t>(jbyt(ites_hdr, 1, 17)) + newpnt;
            const auto nh = static_cast<std::int32_t>(jbyt(ites_hdr, 18, 10));
            const auto idrop = static_cast<std::int32_t>(jbyt(ites_hdr, 28, 5));
            if (idrop != 0 || nh == 0) {
              continue;
            }
            const auto islot = static_cast<std::int32_t>(jbyt(store.iq(loc + iwd1 + 3), 17, 5));
            const std::int32_t locw = loc + iwd1 + 5;
            for (std::int32_t kl = 1; kl <= nh; ++kl) {
              const std::int32_t ites = store.iq(locw + kl);
              const auto nchan = static_cast<std::int32_t>(jbyt(ites, 15, 6));
              if (islot == isldb && nchan == ichdb) {
                timing = static_cast<float>(jbyt(ites, 1, 14));
                found = true;
                break;
              }
            }
          }
        }
        loc = store.lq(loc);
      }
    }
    te.set_float(16, timing);

    // Covariance matrix & timing error
    te.set_float(17, (hit.dx() * hit.dx()) / 3.0f);
    te.set_float(18, 0.0f);
    te.set_float(19, (hit.dy() * hit.dy()) / 3.0f);
    te.set_float(20, 1.0f);

    // Mass assignment
    te.set_float(21, 1.0f);
    te.set_float(22, -1.0f);
    te.set_float(23, 1006.0f);
    te.set_float(24, 0.9f);

    tes.push_back(te);
  }

  return tes;
}

std::int32_t identify_bunch(std::span<const HofTeBank> tes,
                            std::int32_t nbulep,
                            float bunint_ns) {
  constexpr float kTimm = 307.0f;
  constexpr float kWindow = 24.0f;
  constexpr float kFreLtd = 4.0f * 0.117f;
  constexpr std::int32_t kMaxBun = 4;

  const std::int32_t n_bunches = std::clamp(nbulep, 0, kMaxBun);
  std::array<std::int32_t, kMaxBun> nha{};

  for (const auto& te : tes) {
    const float t = te.get_float(16);
    for (std::int32_t ib = 1; ib <= n_bunches; ++ib) {
      const float timexp = kTimm + static_cast<float>(ib - 1) * bunint_ns * kFreLtd;
      if (std::fabs(timexp - t) <= kWindow) {
        nha[static_cast<std::size_t>(ib - 1)] += 1;
      }
    }
  }

  std::int32_t nn0 = 0;
  std::int32_t nb = 0;
  for (std::int32_t ib = 1; ib <= n_bunches; ++ib) {
    if (nha[static_cast<std::size_t>(ib - 1)] > 0) {
      ++nn0;
      nb = ib;
    }
  }
  return (nn0 == 1) ? nb : 0;
}

}  // namespace delphipp::reco::hof
