// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/dst/event_assoc.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace delphipp::dst {

namespace {

float angular_sep(const DstParticle& a, const DstParticle& b) noexcept {
  const float pa = std::sqrt(a.px * a.px + a.py * a.py + a.pz * a.pz);
  const float pb = std::sqrt(b.px * b.px + b.py * b.py + b.pz * b.pz);
  if (pa <= 1e-6f || pb <= 1e-6f) {
    return std::numbers::pi_v<float>;
  }
  const float dot = a.px * b.px + a.py * b.py + a.pz * b.pz;
  return std::acos(std::clamp(dot / (pa * pb), -1.0f, 1.0f));
}

}  // namespace

std::vector<TrackMatch> associate_particles(std::span<const DstParticle> reco,
                                            std::span<const DstParticle> truth,
                                            float max_delta_r) {
  std::vector<TrackMatch> candidates;
  for (std::size_t r = 0; r < reco.size(); ++r) {
    const float pr = std::sqrt(reco[r].px * reco[r].px + reco[r].py * reco[r].py +
                               reco[r].pz * reco[r].pz);
    for (std::size_t t = 0; t < truth.size(); ++t) {
      if (reco[r].charge != truth[t].charge) {
        continue;
      }
      const float dr = angular_sep(reco[r], truth[t]);
      if (dr > max_delta_r) {
        continue;
      }
      const float pt = std::sqrt(truth[t].px * truth[t].px + truth[t].py * truth[t].py +
                                 truth[t].pz * truth[t].pz);
      const float rel_dp = (pt > 1e-4f) ? std::fabs(pr - pt) / pt : 0.0f;
      candidates.push_back(TrackMatch{
          .reco_idx = r,
          .truth_idx = t,
          .delta_r = dr,
          .rel_dp = rel_dp,
      });
    }
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const TrackMatch& a, const TrackMatch& b) {
              if (std::fabs(a.delta_r - b.delta_r) > 1e-4f) {
                return a.delta_r < b.delta_r;
              }
              return a.rel_dp < b.rel_dp;
            });

  std::vector<bool> used_r(reco.size(), false);
  std::vector<bool> used_t(truth.size(), false);
  std::vector<TrackMatch> matched;
  for (const auto& c : candidates) {
    if (!used_r[c.reco_idx] && !used_t[c.truth_idx]) {
      used_r[c.reco_idx] = true;
      used_t[c.truth_idx] = true;
      matched.push_back(c);
    }
  }
  return matched;
}

}  // namespace delphipp::dst
