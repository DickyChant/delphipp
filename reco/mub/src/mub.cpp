// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/mub.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace delphipp::reco::mub {

namespace {

constexpr std::array<std::int32_t, 7> kMaxPos = {5, 4, 5, 4, 3, 4, 3};
constexpr std::array<std::int32_t, 7> kLayP = {0, 5, 9, 14, 18, 21, 25};
constexpr std::array<std::int32_t, 7> kLayPL = {0, 0, 0, 0, 4, 7, 11};

}  // namespace

std::optional<std::int32_t> chamber_index(std::int32_t sector,
                                          std::int32_t layer,
                                          std::int32_t pos) noexcept {
  if (sector < 1 || sector > 50 || layer < 1 || layer > 7 || pos < 1) {
    return std::nullopt;
  }
  const auto lidx = static_cast<std::size_t>(layer - 1);
  if (pos > kMaxPos[lidx] || (sector > 48 && layer <= 3)) {
    return std::nullopt;
  }
  if (sector <= 48) {
    const std::int32_t jcham = kLayP[lidx] + pos;
    return jcham + 28 * (sector - 1);
  }
  const std::int32_t jcham = kLayPL[lidx] + pos;
  return jcham + 1344 + 14 * (sector - 49);
}

std::optional<ChamberLocation> chamber_location(std::int32_t index_1based) noexcept {
  if (index_1based < 1 || index_1based > kNumChambers) {
    return std::nullopt;
  }
  if (index_1based <= 1344) {
    const std::int32_t sector = (index_1based - 1) / 28 + 1;
    const std::int32_t jcham = (index_1based - 1) % 28 + 1;
    for (std::int32_t layer = 7; layer >= 1; --layer) {
      const auto lidx = static_cast<std::size_t>(layer - 1);
      if (jcham > kLayP[lidx]) {
        return ChamberLocation{sector, layer, jcham - kLayP[lidx]};
      }
    }
  } else {
    const std::int32_t rem = index_1based - 1344;
    const std::int32_t sector = (rem - 1) / 14 + 49;
    const std::int32_t jcham = (rem - 1) % 14 + 1;
    for (std::int32_t layer = 7; layer >= 4; --layer) {
      const auto lidx = static_cast<std::size_t>(layer - 1);
      if (jcham > kLayPL[lidx]) {
        return ChamberLocation{sector, layer, jcham - kLayPL[lidx]};
      }
    }
  }
  return std::nullopt;
}

std::int32_t chamber_id(std::int32_t imub,
                        std::int32_t isect,
                        std::int32_t ihemb,
                        std::int32_t icham) noexcept {
  std::int32_t jsect = 0;
  if (ihemb <= 2) {
    jsect = (ihemb - 1) * 24 + isect;
  } else {
    jsect = (isect == 10) ? (47 + ihemb) : (46 + ihemb);
  }
  const std::int32_t jmub = (imub == 1) ? imub : (imub + 1);
  if (jsect <= 48) {
    return icham + 7 * (jmub - 1) + 28 * (jsect - 1);
  }
  return icham + 7 * (imub - 2) + 28 * 48 + 14 * (jsect - 49);
}

float ltd_time_ns(std::int32_t raw_ltd_word,
                  const LtdCalibration& calib) noexcept {
  const auto u = static_cast<std::uint32_t>(raw_ltd_word);
  const std::uint32_t ib0 = u & 1u;
  const std::uint32_t ib1 = (u >> 1) & 1u;

  float tfine = 0.0f;
  if (ib0 == 1u && ib1 == 0u) {
    tfine = 2.0f * calib.fine_step_ns[0];
  } else if (ib0 == 0u && ib1 == 1u) {
    tfine = 4.0f * calib.fine_step_ns[1];
  } else if (ib0 == 1u && ib1 == 1u) {
    tfine = 6.0f * calib.fine_step_ns[2];
  }

  const std::uint32_t jtime = (u & 0x3FFFu) - (2u * ib1) - ib0;
  const float tcors = static_cast<float>(jtime) * calib.coarse_step_ns;
  return tcors + tfine;
}

std::optional<MubTripletResult> reconstruct_triplet(
    std::int32_t layer_1based,
    std::array<float, 3> times_ns,
    const MubChamberCalib& chamber,
    const MubConfig& cfg) noexcept {
  if (layer_1based < 1 || layer_1based > 7) {
    return std::nullopt;
  }
  constexpr float kEpsZ = 10.0f;

  if (cfg.cosmic_mode) {
    for (float& t : times_ns) {
      t -= cfg.cosmic_bcross_ns;
    }
  }

  const float zlen = chamber.length_cm;
  const float zlen2 = 0.5f * zlen;
  const float tof = cfg.tof_ns[static_cast<std::size_t>(layer_1based - 1)];
  const float ct1 = chamber.tzero_ns[0] + tof;
  const float ct2 = chamber.tzero_ns[1] + tof;
  const float ct3 = chamber.tzero_ns[2] + tof;

  const float margin = cfg.is_simulation ? 0.0f : 200.0f;
  if (times_ns[0] < ct1 - margin || times_ns[1] < ct2 - margin ||
      times_ns[2] < ct3 - margin) {
    return std::nullopt;
  }
  if (times_ns[1] < times_ns[0] || times_ns[2] < times_ns[0]) {
    return std::nullopt;
  }
  if (times_ns[0] > cfg.tmax_anode_ns || times_ns[1] > cfg.tmax_dl_ns ||
      times_ns[2] > cfg.tmax_dl_ns) {
    return std::nullopt;
  }

  const float tdif = 0.5f * (times_ns[1] - times_ns[2] - (ct2 - ct3));
  float dl_vel = chamber.delay_poly[1];
  if (cfg.is_simulation && dl_vel <= 0.0f) {
    dl_vel = 0.547f;
  }

  float zloc = 0.0f;
  if (cfg.is_simulation) {
    zloc = dl_vel * tdif;
  } else {
    const float td2 = tdif * tdif;
    zloc = chamber.delay_poly[0] + dl_vel * tdif + chamber.delay_poly[2] * td2 +
           chamber.delay_poly[3] * td2 * tdif + chamber.delay_poly[4] * td2 * td2;
  }

  if (std::fabs(zloc) > zlen2 + kEpsZ) {
    return std::nullopt;
  }
  zloc = std::clamp(zloc, -zlen2, zlen2);

  float tpr = (zloc + zlen2) / cfg.vanode_cm_per_ns;
  if (tpr < 0.0f) {
    tpr = 0.0f;
  }

  const float tdd = times_ns[0] - ct1 - tpr;
  const float td = (tdd <= 0.0f) ? 0.0f : tdd;

  float ddrft = 0.0f;
  if (cfg.is_simulation) {
    ddrft = chamber.anode_poly[1] * td;
  } else {
    float velcor = cfg.veldrf_cm_per_ns;
    if (cfg.atmosp_mbar > 960.0f && cfg.atmosp_mbar < 1000.0f) {
      velcor = cfg.veldrf_cm_per_ns * 980.0f / cfg.atmosp_mbar;
    }
    ddrft = (td < 50.0f) ? (1.5f * velcor * td) : (0.117f + velcor * td);
  }
  ddrft = std::clamp(ddrft, 0.0f, chamber.half_width_cm);

  const float tsum =
      times_ns[1] + times_ns[2] - 2.0f * times_ns[0] - (ct2 + ct3 - 2.0f * ct1) + 2.0f * tpr;

  return MubTripletResult{
      .z_local_cm = -zloc + zlen2,
      .drift_cm = ddrft,
      .anode_time_ns = tdd,
      .time_sum_ns = tsum,
  };
}

std::vector<MubRawHit> unpack_raw_blocklet(std::span<const std::int32_t> words,
                                           const LtdCalibration& ltd) {
  std::vector<MubRawHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  LtdCalibration eff = ltd;
  if (eff.coarse_step_ns == 0.0f) {
    eff.coarse_step_ns = 8.0f;
  }
  constexpr float kTwoPi = 2.0f * std::numbers::pi_v<float>;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    if ((u & 0x3FFFu) == 0u) {
      continue;
    }
    const std::int32_t sector = static_cast<std::int32_t>((u >> 22) & 0x3Fu) % 48 + 1;
    const std::int32_t layer = static_cast<std::int32_t>((u >> 19) & 0x7u) % 7 + 1;
    const std::int32_t pos = static_cast<std::int32_t>((u >> 16) & 0x7u) % 3 + 1;
    const float r_cm = 430.0f + static_cast<float>(layer) * 10.0f;
    const float phi_rad = (static_cast<float>((sector - 1) % 24) + 0.5f) * (kTwoPi / 24.0f);
    const float z_cm = (sector <= 24) ? 100.0f : -100.0f;
    hits.push_back(MubRawHit{
        .sector = sector,
        .layer = layer,
        .pos = pos,
        .r_cm = r_cm,
        .phi_rad = phi_rad,
        .z_cm = z_cm,
        .time_ns = ltd_time_ns(words[i], eff),
    });
  }
  return hits;
}

std::vector<TanagraId> write_ter_banks(TanagraStore& store,
                                       std::span<const MubRawHit> hits) {
  std::vector<TanagraId> ids;
  ids.reserve(hits.size());
  for (const auto& h : hits) {
    ids.push_back(store.add_ter(TerBank{
        .detector_id = kDetectorId,
        .submod_id = h.sector * 10 + h.layer,
        .descriptor = 125,
        .fwords = {h.r_cm, h.r_cm * h.phi_rad, h.z_cm, h.time_ns},
    }));
  }
  return ids;
}

}  // namespace delphipp::reco::mub
