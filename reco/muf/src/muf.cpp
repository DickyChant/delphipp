// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/muf.hpp"

#include <cmath>

namespace delphipp::reco::muf {

float time_slew(int end_1or2,
                float chamber_len_cm,
                float vline_cm_per_ns,
                float time_ns,
                const MufConfig& cfg) noexcept {
  if (cfg.is_simulation) {
    return 0.0f;
  }
  const std::size_t k = (chamber_len_cm < 400.0f) ? 1 : 0;
  if (end_1or2 == 1) {
    const float x = (time_ns * vline_cm_per_ns - chamber_len_cm * 0.5f) / 100.0f;
    return cfg.ctslw0[k] + cfg.ctslw1[k] * x + cfg.ctslw2[k] * x * x +
           cfg.ctslw3[k] * x * x * x;
  }
  if (end_1or2 == 2) {
    const float x = (chamber_len_cm * 0.5f - time_ns * vline_cm_per_ns) / 100.0f;
    return cfg.ctslw0[k] - cfg.ctslw1[k] * x + cfg.ctslw2[k] * x * x -
           cfg.ctslw3[k] * x * x * x;
  }
  return 0.0f;
}

DelayLineResult delay_line_distance(const MufRawTriplet& h,
                                    float vanode,
                                    const MufConfig& cfg) noexcept {
  DelayLineResult res{};
  const float vl = h.vline_cm_per_ns;
  const float chl = h.chamber_len_cm;

  auto end1_only = [&]() {
    const float v = vanode * vl / (vanode - vl);
    res.x_on_line_cm = v * (h.tl1_ns - h.ta_ns) - v * (h.tc1_ns - h.tca_ns);
    res.sigma_cm = 1.414f * cfg.sigma_delay_line_cm;
  };
  auto end2_only = [&]() {
    const float v = vanode * vl / (vanode + vl);
    res.x_on_line_cm =
        chl * vanode / (vl + vanode) + v * (h.tc2_ns - h.tca_ns) - v * (h.tl2_ns - h.ta_ns);
    res.sigma_cm = 1.414f * cfg.sigma_delay_line_cm;
  };

  if (h.tl1_ns > 0.0f && h.tl2_ns > 0.0f) {
    res.x_on_line_cm =
        chl * 0.5f - vl * (h.tc1_ns - h.tc2_ns) * 0.5f + vl * (h.tl1_ns - h.tl2_ns) * 0.5f;
    res.sigma_cm = cfg.sigma_delay_line_cm;
    if (h.time_sum_ns <= cfg.prepulse_cut_ns) {
      if (res.x_on_line_cm >= chl * 0.5f) {
        end1_only();
      } else {
        end2_only();
      }
    }
  } else if (h.tl1_ns > 0.0f) {
    end1_only();
  } else if (h.tl2_ns > 0.0f) {
    end2_only();
  }

  return res;
}

float wire_position(std::int32_t quadrant_1based,
                    std::int32_t chamber_1based,
                    std::int32_t layer_1or2,
                    float x_on_line_cm,
                    const MufCalibration& calib) noexcept {
  const auto q = static_cast<std::size_t>(quadrant_1based - 1);
  const float step = static_cast<float>(chamber_1based - 1) * calib.wire_pitch_cm;
  if (layer_1or2 == 1) {
    const float angle = calib.qgeom[q][4];
    const float axis = calib.qgeom[q][2] + step / std::cos(angle);
    return axis - x_on_line_cm * std::sin(angle);
  }
  const float angle = calib.qgeom[q][5];
  const float axis = calib.qgeom[q][3] + step / std::cos(angle);
  return axis - x_on_line_cm * std::sin(angle);
}

float drift_distance(std::int32_t quadrant_1based,
                     std::int32_t layer_1or2,
                     float drift_time_ns,
                     const MufCalibration& calib) noexcept {
  const auto q = static_cast<std::size_t>(quadrant_1based - 1);
  const auto l = static_cast<std::size_t>(layer_1or2 - 1);
  const float x0 = calib.sx0par[q][l];
  const float vdrift = calib.svdrif[q][l];
  const float t1 = x0 / (2.0f * vdrift);
  if (drift_time_ns < t1) {
    return x0 * std::sqrt(drift_time_ns / t1);
  }
  return x0 + vdrift * (drift_time_ns - t1);
}

std::vector<MufSpaceHit> reconstruct_hits(std::span<const MufRawTriplet> triplets,
                                          const MufCalibration& calib,
                                          const MufConfig& cfg) {
  std::vector<MufSpaceHit> hits;
  hits.reserve(triplets.size());

  for (const auto& t : triplets) {
    if (t.quadrant < 1 || t.quadrant > static_cast<std::int32_t>(kNumQuadrants) ||
        (t.layer != 1 && t.layer != 2)) {
      continue;
    }

    const auto dl = delay_line_distance(t, calib.vanode_cm_per_ns, cfg);
    const float xwire = wire_position(t.quadrant, t.chamber, t.layer, dl.x_on_line_cm, calib);

    const float tpro = dl.x_on_line_cm / calib.vanode_cm_per_ns;
    float tdr = t.ta_ns - cfg.flight_time_ns - tpro - t.tca_ns;
    if (tdr <= 0.0f) {
      tdr = 1e-6f;
    }
    const float dtow = drift_distance(t.quadrant, t.layer, tdr, calib);

    hits.push_back(MufSpaceHit{
        .quadrant = t.quadrant,
        .chamber = t.chamber,
        .layer = t.layer,
        .encoded_id = t.quadrant * 1000 + t.chamber * 10 + t.layer,
        .drift_dist_cm = dtow,
        .x_on_line_cm = dl.x_on_line_cm,
        .x_wire_cm = xwire,
        .sigma_drift_cm = cfg.sigma_drift_cm,
        .sigma_line_cm = dl.sigma_cm,
    });
  }

  return hits;
}

}  // namespace delphipp::reco::muf
