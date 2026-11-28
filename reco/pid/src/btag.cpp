// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/pid/btag.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "pid_stats.hpp"

namespace delphipp::reco::pid {

namespace {

[[nodiscard]] Vec3 normalize_vec(const Vec3& v, const Vec3& fallback) noexcept {
  const double n2 = static_cast<double>(v[0]) * v[0] +
                    static_cast<double>(v[1]) * v[1] +
                    static_cast<double>(v[2]) * v[2];
  if (n2 > 1e-12) {
    const double inv = 1.0 / std::sqrt(n2);
    return Vec3{
        static_cast<float>(v[0] * inv),
        static_cast<float>(v[1] * inv),
        static_cast<float>(v[2] * inv),
    };
  }
  const double f2 = static_cast<double>(fallback[0]) * fallback[0] +
                    static_cast<double>(fallback[1]) * fallback[1] +
                    static_cast<double>(fallback[2]) * fallback[2];
  if (f2 > 1e-12) {
    const double inv = 1.0 / std::sqrt(f2);
    return Vec3{
        static_cast<float>(fallback[0] * inv),
        static_cast<float>(fallback[1] * inv),
        static_cast<float>(fallback[2] * inv),
    };
  }
  return Vec3{0.0f, 0.0f, 1.0f};
}

[[nodiscard]] double dot3(const Vec3& a, const Vec3& b) noexcept {
  return static_cast<double>(a[0]) * b[0] +
         static_cast<double>(a[1]) * b[1] +
         static_cast<double>(a[2]) * b[2];
}

}  // namespace

float track_resolution_probability(float abs_significance,
                                   std::int32_t n_vd_layers,
                                   std::span<const float> custom_table_400) noexcept {
  const double s = std::max(0.0, static_cast<double>(std::fabs(abs_significance)));

  // Optional AATRKP 400-bin table with quadratic interpolation (aabtagxx.car lines 3477-3500)
  if (custom_table_400.size() >= 3) {
    constexpr double kMaxS = 40.0;
    const double step = kMaxS / static_cast<double>(custom_table_400.size());
    const std::size_t n = custom_table_400.size();
    const std::size_t ib = static_cast<std::size_t>(s / step);
    if (ib == 0) {
      const double p1 = custom_table_400[0];
      const double p2 = custom_table_400[1];
      return static_cast<float>(std::clamp(p1 + (p2 - p1) * (s / step), 1e-12, 1.0));
    }
    if (ib + 1 < n) {
      const double p1 = custom_table_400[ib];
      const double p2 = custom_table_400[ib + 1];
      const double p3 = custom_table_400[ib - 1];
      const double x0 = static_cast<double>(ib) * step;
      const double a = 0.5 * (p2 + p3 - 2.0 * p1) / (step * step);
      const double b = (p2 - p3 - 4.0 * a * x0 * step) / (2.0 * step);
      const double c = p1 - a * x0 * x0 - b * x0;
      return static_cast<float>(std::clamp(a * s * s + b * s + c, 1e-12, 1.0));
    }
    return static_cast<float>(std::clamp(static_cast<double>(custom_table_400.back()), 1e-12, 1.0));
  }

  // Calibrated DELPHI 94D VD resolution function P_T(S) = integral_{|S|}^inf R(x) dx:
  // Gaussian core + modest non-Gaussian tracking/scattering tails (DELPHI 94-125).
  double w1 = 0.92;
  double s1 = 1.00;
  double w2 = 0.065;
  double s2 = 1.75;
  double w3 = 0.015;
  double s3 = 3.20;

  if (n_vd_layers == 2) {
    w1 = 0.89; s1 = 1.04;
    w2 = 0.09; s2 = 1.85;
    w3 = 0.02; s3 = 3.40;
  } else if (n_vd_layers <= 1) {
    w1 = 0.85; s1 = 1.10;
    w2 = 0.12; s2 = 2.00;
    w3 = 0.03; s3 = 3.60;
  }

  const double p =
      w1 * std::erfc(s / (std::numbers::sqrt2 * s1)) +
      w2 * std::erfc(s / (std::numbers::sqrt2 * s2)) +
      w3 * std::erfc(s / (std::numbers::sqrt2 * s3));

  return static_cast<float>(std::clamp(p, 1e-12, 1.0));
}

float fisher_combine_probabilities(std::span<const float> track_probs) noexcept {
  std::vector<double> vals;
  vals.reserve(track_probs.size());
  for (float p : track_probs) {
    if (p > 0.0f && p < 1.0f) {
      vals.push_back(static_cast<double>(p));
    }
  }
  return static_cast<float>(detail::fisher_combined_prob(vals));
}

BTagTrackResult evaluate_track_btag(
    const BTagTrackInput& track_in,
    const reco::tracking::VertexFitResult& primary_vertex,
    const Vec3& thrust_axis,
    const BTagConfig& cfg) noexcept {
  BTagTrackResult res{};

  const PerigeeTrack& trk = track_in.track;
  const float p_gev = trk.p_gev(cfg.bz_gev_per_cm);
  if (!std::isfinite(p_gev) || p_gev < cfg.min_momentum_gev ||
      track_in.n_vd_rphi < 1) {
    return res;
  }

  const Vec3& pv = primary_vertex.xyz_cm;
  const double sin_th = std::max(std::fabs(std::sin(static_cast<double>(trk.theta))), 1e-4);
  const double cos_th = std::cos(static_cast<double>(trk.theta));
  const double cos_phi = std::cos(static_cast<double>(trk.phi0));
  const double sin_phi = std::sin(static_cast<double>(trk.phi0));

  const Vec3 dir_trk{
      static_cast<float>(sin_th * cos_phi),
      static_cast<float>(sin_th * sin_phi),
      static_cast<float>(cos_th),
  };

  const Vec3 thrust = normalize_vec(thrust_axis, Vec3{0.0f, 0.0f, 1.0f});
  const Vec3 jet = normalize_vec(track_in.jet_axis, thrust);
  res.hemisphere = (dot3(dir_trk, thrust) >= 0.0) ? 1 : 2;

  // Impact parameter wrt primary vertex in r-phi and z
  // Perigee point: x0 = -eps*sin(phi0), y0 = eps*cos(phi0), z0 = trk.z0
  const double x0 = -static_cast<double>(trk.eps) * sin_phi;
  const double y0 = static_cast<double>(trk.eps) * cos_phi;
  const double z0 = static_cast<double>(trk.z0);

  // Linearized s_xy to closest approach to PV
  const double s_xy_pv =
      (static_cast<double>(pv[0]) - x0) * cos_phi +
      (static_cast<double>(pv[1]) - y0) * sin_phi;
  const auto pt_pv = propagate_arc(trk, static_cast<float>(s_xy_pv));

  const double dx_pv = static_cast<double>(pt_pv.xyz[0]) - pv[0];
  const double dy_pv = static_cast<double>(pt_pv.xyz[1]) - pv[1];
  const double dz_pv = static_cast<double>(pt_pv.xyz[2]) - pv[2];

  // Signed r-phi distance along n_rphi = (-sin(phi0), cos(phi0))
  const double d_rphi_raw = -dx_pv * sin_phi + dy_pv * cos_phi;
  const double d_z_perp = dz_pv * sin_th;
  const double d_3d_mag = std::hypot(d_rphi_raw, d_z_perp);

  // Closest approach between track line/helix and jet axis ray pv + u * jet (AADIST)
  const Vec3 w0{
      static_cast<float>(x0 - pv[0]),
      static_cast<float>(y0 - pv[1]),
      static_cast<float>(z0 - pv[2]),
  };
  const double b = dot3(dir_trk, jet);
  const double denom = 1.0 - b * b;

  double u_jet = 0.0;
  double s3d_jet = 0.0;
  if (denom > 1e-6) {
    const double w0_j = dot3(w0, jet);
    const double w0_t = dot3(w0, dir_trk);
    u_jet = (w0_j - b * w0_t) / denom;
    s3d_jet = b * u_jet - w0_t;
  } else {
    u_jet = -dot3(w0, jet);
  }

  // Refine 1 step along exact helix around s_xy = s3d_jet * sin_th
  const auto pt_jet = propagate_arc(trk, static_cast<float>(s3d_jet * sin_th));
  const Vec3 pt_rel{
      pt_jet.xyz[0] - pv[0],
      pt_jet.xyz[1] - pv[1],
      pt_jet.xyz[2] - pv[2],
  };
  const Vec3 pt_dir = pt_jet.dir;
  const double b_ref = dot3(pt_dir, jet);
  const double denom_ref = 1.0 - b_ref * b_ref;
  if (denom_ref > 1e-6) {
    u_jet = (dot3(pt_rel, jet) - b_ref * dot3(pt_rel, pt_dir)) / denom_ref;
  }

  const Vec3 jet_closest{
      pv[0] + static_cast<float>(u_jet) * jet[0],
      pv[1] + static_cast<float>(u_jet) * jet[1],
      pv[2] + static_cast<float>(u_jet) * jet[2],
  };
  const double dist_to_jet = std::sqrt(
      std::pow(static_cast<double>(pt_jet.xyz[0]) - jet_closest[0], 2.0) +
      std::pow(static_cast<double>(pt_jet.xyz[1]) - jet_closest[1], 2.0) +
      std::pow(static_cast<double>(pt_jet.xyz[2]) - jet_closest[2], 2.0));
  res.dist_to_jet_cm = static_cast<float>(dist_to_jet);

  // AADCAJ sign convention: +1 if DCA to jet axis lies downstream (u_jet >= 0), else -1
  const float sign = (u_jet >= 0.0) ? 1.0f : -1.0f;

  // Combine track perigee covariance and primary-vertex uncertainty in r-phi and z
  const double pv_var_rphi =
      sin_phi * sin_phi * static_cast<double>(primary_vertex.sigma_cm[0]) *
          static_cast<double>(primary_vertex.sigma_cm[0]) +
      cos_phi * cos_phi * static_cast<double>(primary_vertex.sigma_cm[1]) *
          static_cast<double>(primary_vertex.sigma_cm[1]);
  const double pv_var_z =
      static_cast<double>(primary_vertex.sigma_cm[2]) *
      static_cast<double>(primary_vertex.sigma_cm[2]);

  const double sig_rphi = std::sqrt(
      std::pow(static_cast<double>(track_in.sigma_rphi_cm), 2.0) +
      std::max(0.0, pv_var_rphi));
  const double sig_z = std::sqrt(
      std::pow(static_cast<double>(track_in.sigma_z_cm), 2.0) +
      std::max(0.0, pv_var_z));

  res.sigma_rphi_cm = static_cast<float>(std::max(sig_rphi, static_cast<double>(cfg.min_sigma_cm)));
  res.sigma_z_cm = static_cast<float>(std::max(sig_z, static_cast<double>(cfg.min_sigma_cm)));

  const double w_rphi = std::fabs(d_rphi_raw);
  const double w_z = std::fabs(d_z_perp);
  const double sig_z_perp = res.sigma_z_cm * sin_th;
  double sig_3d = res.sigma_rphi_cm;
  if (d_3d_mag > 1e-9 && cfg.use_3d_z && track_in.n_vd_z > 0) {
    sig_3d = std::sqrt(
        (w_rphi * w_rphi * res.sigma_rphi_cm * res.sigma_rphi_cm +
         w_z * w_z * sig_z_perp * sig_z_perp) /
        (d_3d_mag * d_3d_mag));
  }
  res.sigma_3d_cm = static_cast<float>(std::max(sig_3d, static_cast<double>(cfg.min_sigma_cm)));

  res.ip_rphi_cm = sign * static_cast<float>(std::fabs(d_rphi_raw));
  res.ip_z_cm = sign * static_cast<float>(std::fabs(dz_pv));
  res.ip_3d_cm = sign * static_cast<float>(d_3d_mag);

  res.significance_rphi = res.ip_rphi_cm / res.sigma_rphi_cm;
  res.significance_z = res.ip_z_cm / res.sigma_z_cm;
  res.significance_3d = res.ip_3d_cm / res.sigma_3d_cm;

  // Selection cuts (AARJCT / AATPRB)
  if (std::fabs(d_rphi_raw) > cfg.max_impact_rphi_cm) {
    return res;
  }
  if (cfg.use_3d_z && track_in.n_vd_z > 0 && std::fabs(dz_pv) > cfg.max_impact_z_cm) {
    return res;
  }
  if (dist_to_jet / res.sigma_rphi_cm > cfg.max_dist_to_jet_sig &&
      dist_to_jet > 0.25) {
    return res;
  }

  res.prob_rphi = track_resolution_probability(res.significance_rphi, track_in.n_vd_rphi);
  if (cfg.use_3d_z && track_in.n_vd_z > 0) {
    res.prob_z = track_resolution_probability(res.significance_z, track_in.n_vd_z);
    res.prob_3d = track_resolution_probability(res.significance_3d, track_in.n_vd_rphi);
  } else {
    res.prob_z = 1.0f;
    res.prob_3d = res.prob_rphi;
  }
  res.selected = true;
  return res;
}

BTagEventResult evaluate_event_btag(
    std::span<const BTagTrackInput> tracks,
    const reco::tracking::VertexFitResult& primary_vertex,
    const Vec3& thrust_axis,
    const BTagConfig& cfg) noexcept {
  BTagEventResult out{};
  out.tracks.reserve(tracks.size());

  std::vector<double> pos_probs;
  std::vector<double> neg_probs;
  std::vector<double> all_probs;
  std::array<std::vector<double>, 2> hemi_pos;
  std::array<std::vector<double>, 2> hemi_neg;

  for (const auto& t : tracks) {
    const auto tr = evaluate_track_btag(t, primary_vertex, thrust_axis, cfg);
    out.tracks.push_back(tr);
    if (!tr.selected) {
      continue;
    }
    ++out.num_selected_tracks;

    const std::size_t hi = (tr.hemisphere == 1) ? 0u : 1u;
    const bool is_pos = (tr.significance_3d >= 0.0f);
    if (is_pos) {
      ++out.num_positive_tracks;
    } else {
      ++out.num_negative_tracks;
    }

    auto push_prob = [&](float p) {
      if (p > 0.0f && p < 1.0f) {
        const double dp = static_cast<double>(p);
        all_probs.push_back(dp);
        if (is_pos) {
          pos_probs.push_back(dp);
          hemi_pos[hi].push_back(dp);
        } else {
          neg_probs.push_back(dp);
          hemi_neg[hi].push_back(dp);
        }
      }
    };

    // Combine r-phi (and z when VD z layers present), matching bbPROB lines 3620-3634
    push_prob(tr.prob_rphi);
    if (cfg.use_3d_z && t.n_vd_z > 0) {
      push_prob(tr.prob_z);
    }
  }

  out.p_event_btag = static_cast<float>(detail::fisher_combined_prob(pos_probs));
  out.p_event_negative = static_cast<float>(detail::fisher_combined_prob(neg_probs));
  out.p_event_all = static_cast<float>(detail::fisher_combined_prob(all_probs));
  out.p_hemisphere_btag[0] = static_cast<float>(detail::fisher_combined_prob(hemi_pos[0]));
  out.p_hemisphere_btag[1] = static_cast<float>(detail::fisher_combined_prob(hemi_pos[1]));
  out.p_hemisphere_negative[0] = static_cast<float>(detail::fisher_combined_prob(hemi_neg[0]));
  out.p_hemisphere_negative[1] = static_cast<float>(detail::fisher_combined_prob(hemi_neg[1]));

  return out;
}

}  // namespace delphipp::reco::pid
