// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/tof.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace delphipp::reco::tof {

namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kTwoPi = 2.0f * kPi;
constexpr std::int32_t kIgDesc = 125;  // cylindrical (R, Rphi, z)
constexpr float kPtDesc = 14.0f;

}  // namespace

TofSectorPose sector_pose(std::int32_t sector_1based) noexcept {
  const float z = (sector_1based > 24) ? -11.5f : 11.5f;
  float phi = (-15.0f * static_cast<float>(sector_1based - 1) - 7.5f + 90.0f) *
              (kPi / 180.0f);
  while (phi < 0.0f) {
    phi += kTwoPi;
  }
  while (phi >= kTwoPi) {
    phi -= kTwoPi;
  }
  return TofSectorPose{.r_cm = kTofRadiusCm, .phi_rad = phi, .z_cm = z};
}

std::vector<TofCounterHit> unpack_raw_blocklet(std::span<const std::int32_t> words) {
  std::vector<TofCounterHit> hits;
  if (words.size() < 4) {
    return hits;
  }
  const std::size_t start =
      (words[0] >= 4 && static_cast<std::size_t>(words[0]) <= words.size()) ? 4u : 0u;
  for (std::size_t i = start; i < words.size(); ++i) {
    const std::uint32_t u = static_cast<std::uint32_t>(words[i]);
    const std::int32_t slot = static_cast<std::int32_t>((u >> 27) & 0x1Fu);
    const std::int32_t chan = static_cast<std::int32_t>((u >> 17) & 0x7Fu);
    const std::uint32_t tdc = u & 0xFFFFu;
    if (slot <= 0 || chan >= 64 || tdc == 0u || tdc >= 0xFFFFu) {
      continue;
    }
    const std::int32_t counter_num = ((slot * 16 + chan) % 192) + 1;
    const std::int32_t sector = ((counter_num - 1) / 4) + 1;
    const auto pose = sector_pose(sector);
    const float time_ns = 0.5f * static_cast<float>(tdc);
    hits.push_back(TofCounterHit{
        .counter_num = counter_num,
        .r_cm = pose.r_cm,
        .r_phi_cm = pose.r_cm * pose.phi_rad,
        .z_cm = pose.z_cm,
        .pm_ext = time_ns,
        .pm_int = time_ns,
    });
  }
  return hits;
}

std::vector<TanagraId> create_td_banks(
    TanagraStore& store,
    std::span<const TofCounterHit> hits,
    std::span<const std::int32_t> multihit_sectors,
    std::int32_t year_version) {
  std::array<bool, kNumSectors + 1> seen_sector{};
  std::vector<TanagraId> ids;

  for (const auto& h : hits) {
    const std::int32_t nsect = (h.counter_num + 3) / 4;
    if (nsect >= 1 && nsect <= kNumSectors) {
      seen_sector[static_cast<std::size_t>(nsect)] = true;
    }
    TdrBank td;
    td.detector_id = kDetectorId;
    td.submod_id = h.counter_num;
    td.descriptor = kIgDesc;
    td.num_points = 1;
    td.words_per_point = 9;
    td.header_floats = {1.0f, 10.0f, 35.50f};
    td.data = {kPtDesc, h.r_cm, h.r_phi_cm, h.z_cm, 1.0f, 10.0f, 35.50f, h.pm_ext, h.pm_int};
    ids.push_back(store.add_tdr(std::move(td)));
  }

  if (year_version < 95) {
    for (const std::int32_t ilis : multihit_sectors) {
      if (ilis < 1 || ilis > kNumSectors) {
        continue;
      }
      if (seen_sector[static_cast<std::size_t>(ilis)]) {
        continue;
      }
      seen_sector[static_cast<std::size_t>(ilis)] = true;
      const auto pose = sector_pose(ilis);
      TdrBank td;
      td.detector_id = kDetectorId;
      td.submod_id = 192 + 100 * ilis;
      td.descriptor = kIgDesc;
      td.num_points = 1;
      td.words_per_point = 9;
      td.header_floats = {1.0f, 20.0f, 360.0f};
      td.data = {kPtDesc, pose.r_cm, pose.r_cm * pose.phi_rad, pose.z_cm,
                 1.0f, 20.0f, 360.0f, 0.0f, 0.0f};
      ids.push_back(store.add_tdr(std::move(td)));
    }
  }

  return ids;
}

std::int32_t identify_bunch_majority(std::span<const float> tmaj2_times_ns,
                                     std::int32_t nbulep,
                                     float bunint_ns) noexcept {
  if (tmaj2_times_ns.empty() || nbulep <= 0) {
    return 0;
  }
  float rmin = 2000.0f;
  for (const float t : tmaj2_times_ns) {
    rmin = std::min(rmin, t);
  }
  std::int32_t nbunch = 0;
  constexpr float kXMin1 = 730.0f;
  constexpr float kXMax1 = 800.0f;
  for (std::int32_t i = 1; i <= nbulep; ++i) {
    const float offset = static_cast<float>(i - 1) * bunint_ns;
    if (rmin >= kXMin1 + offset && rmin <= kXMax1 + offset) {
      nbunch = i;
    }
  }
  return nbunch;
}

}  // namespace delphipp::reco::tof
