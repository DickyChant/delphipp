// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/v0_finder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace delphipp::reco::tracking {

namespace {

float inv_mass_two_body(const Vec3& p1, float m1,
                        const Vec3& p2, float m2) noexcept {
  const double p1_2 = static_cast<double>(p1[0]) * p1[0] +
                      static_cast<double>(p1[1]) * p1[1] +
                      static_cast<double>(p1[2]) * p1[2];
  const double p2_2 = static_cast<double>(p2[0]) * p2[0] +
                      static_cast<double>(p2[1]) * p2[1] +
                      static_cast<double>(p2[2]) * p2[2];
  const double e1 = std::sqrt(p1_2 + static_cast<double>(m1) * m1);
  const double e2 = std::sqrt(p2_2 + static_cast<double>(m2) * m2);
  const double px = static_cast<double>(p1[0]) + p2[0];
  const double py = static_cast<double>(p1[1]) + p2[1];
  const double pz = static_cast<double>(p1[2]) + p2[2];
  const double m2_tot = (e1 + e2) * (e1 + e2) - (px * px + py * py + pz * pz);
  return static_cast<float>(std::sqrt(std::max(0.0, m2_tot)));
}

float dist3(const Vec3& a, const Vec3& b) noexcept {
  const float dx = a[0] - b[0];
  const float dy = a[1] - b[1];
  const float dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

}  // namespace

std::vector<V0Candidate> find_v0_candidates(
    std::span<const PerigeeTrack> tracks,
    std::span<const std::int32_t> charges,
    float bz_gev_per_cm,
    const Vec3& pv_xyz_cm,
    float max_dca_cm) {
  std::vector<V0Candidate> out;
  const std::size_t n = std::min(tracks.size(), charges.size());
  if (n < 2) {
    return out;
  }

  for (std::size_t ip = 0; ip < n; ++ip) {
    if (charges[ip] <= 0) {
      continue;
    }
    const float p_pos = tracks[ip].p_gev(bz_gev_per_cm);
    if (!std::isfinite(p_pos) || p_pos <= 0.05f || p_pos > 200.0f) {
      continue;
    }

    for (std::size_t in = 0; in < n; ++in) {
      if (charges[in] >= 0) {
        continue;
      }
      const float p_neg = tracks[in].p_gev(bz_gev_per_cm);
      if (!std::isfinite(p_neg) || p_neg <= 0.05f || p_neg > 200.0f) {
        continue;
      }

      // Coarse + refined 2D search in transverse arc lengths (s_pos, s_neg)
      float best_s1 = 0.0f;
      float best_s2 = 0.0f;
      float best_d = std::numeric_limits<float>::infinity();

      for (int i1 = 0; i1 <= 30; ++i1) {
        const float s1 = static_cast<float>(i1) * 2.0f;
        const auto pt1 = propagate_arc(tracks[ip], s1);
        for (int i2 = 0; i2 <= 30; ++i2) {
          const float s2 = static_cast<float>(i2) * 2.0f;
          const auto pt2 = propagate_arc(tracks[in], s2);
          const float d = dist3(pt1.xyz, pt2.xyz);
          if (d < best_d) {
            best_d = d;
            best_s1 = s1;
            best_s2 = s2;
          }
        }
      }

      for (float step : {0.25f, 0.025f}) {
        const float c1 = best_s1;
        const float c2 = best_s2;
        for (int k1 = -8; k1 <= 8; ++k1) {
          const float s1 = std::max(0.0f, c1 + static_cast<float>(k1) * step);
          const auto pt1 = propagate_arc(tracks[ip], s1);
          for (int k2 = -8; k2 <= 8; ++k2) {
            const float s2 = std::max(0.0f, c2 + static_cast<float>(k2) * step);
            const auto pt2 = propagate_arc(tracks[in], s2);
            const float d = dist3(pt1.xyz, pt2.xyz);
            if (d < best_d) {
              best_d = d;
              best_s1 = s1;
              best_s2 = s2;
            }
          }
        }
      }

      if (best_d > max_dca_cm) {
        continue;
      }

      const auto pt1 = propagate_arc(tracks[ip], best_s1);
      const auto pt2 = propagate_arc(tracks[in], best_s2);
      const Vec3 vtx{
          0.5f * (pt1.xyz[0] + pt2.xyz[0]),
          0.5f * (pt1.xyz[1] + pt2.xyz[1]),
          0.5f * (pt1.xyz[2] + pt2.xyz[2]),
      };

      const Vec3 p1_vec{
          p_pos * pt1.dir[0],
          p_pos * pt1.dir[1],
          p_pos * pt1.dir[2],
      };
      const Vec3 p2_vec{
          p_neg * pt2.dir[0],
          p_neg * pt2.dir[1],
          p_neg * pt2.dir[2],
      };

      const Vec3 pv_to_sv{
          vtx[0] - pv_xyz_cm[0],
          vtx[1] - pv_xyz_cm[1],
          vtx[2] - pv_xyz_cm[2],
      };
      const Vec3 p_tot{
          p1_vec[0] + p2_vec[0],
          p1_vec[1] + p2_vec[1],
          p1_vec[2] + p2_vec[2],
      };
      const float flight_dot = pv_to_sv[0] * p_tot[0] +
                               pv_to_sv[1] * p_tot[1] +
                               pv_to_sv[2] * p_tot[2];
      const float decay_len = std::sqrt(pv_to_sv[0] * pv_to_sv[0] +
                                        pv_to_sv[1] * pv_to_sv[1] +
                                        pv_to_sv[2] * pv_to_sv[2]);
      if (decay_len > 0.2f && flight_dot < 0.0f) {
        continue;
      }

      constexpr float kMe = 0.000511f;
      constexpr float kMpi = 0.139570f;
      constexpr float kMp = 0.938272f;
      const float m_ks = inv_mass_two_body(p1_vec, kMpi, p2_vec, kMpi);
      const float m_lam = inv_mass_two_body(p1_vec, kMp, p2_vec, kMpi);
      const float m_alam = inv_mass_two_body(p1_vec, kMpi, p2_vec, kMp);
      const float m_ee = inv_mass_two_body(p1_vec, kMe, p2_vec, kMe);

      V0Hypothesis hyp = V0Hypothesis::KShort;
      if (m_ee < 0.08f) {
        hyp = V0Hypothesis::PhotonConversion;
      } else {
        const float d_ks = std::fabs(m_ks - 0.497611f);
        const float d_lam = std::fabs(m_lam - 1.115683f);
        const float d_alam = std::fabs(m_alam - 1.115683f);
        if (d_lam <= d_ks && d_lam <= d_alam) {
          hyp = V0Hypothesis::Lambda;
        } else if (d_alam <= d_ks && d_alam < d_lam) {
          hyp = V0Hypothesis::AntiLambda;
        }
      }

      out.push_back(V0Candidate{
          .pos_idx = ip,
          .neg_idx = in,
          .decay_xyz_cm = vtx,
          .decay_length_cm = decay_len,
          .dca_cm = best_d,
          .mass_kshort_gev = m_ks,
          .mass_lambda_gev = m_lam,
          .mass_antilambda_gev = m_alam,
          .mass_ee_gev = m_ee,
          .best_hypothesis = hyp,
      });
    }
  }
  return out;
}

}  // namespace delphipp::reco::tracking
