// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/vertex_fitter.hpp"

#include <cmath>

namespace delphipp::reco::tracking {

VertexFitResult fit_primary_vertex(std::span<const PerigeeTrack> tracks,
                                   const BeamspotConstraint& beamspot,
                                   float track_sigma_rphi_cm,
                                   float track_sigma_z_cm) noexcept {
  // Diagonal-beamspot + linearized perigee normal equations for (x, y) and z
  double a11 = 1.0 / (static_cast<double>(beamspot.sigma_cm[0]) * beamspot.sigma_cm[0]);
  double a22 = 1.0 / (static_cast<double>(beamspot.sigma_cm[1]) * beamspot.sigma_cm[1]);
  double a12 = 0.0;
  double bx = a11 * static_cast<double>(beamspot.xyz_cm[0]);
  double by = a22 * static_cast<double>(beamspot.xyz_cm[1]);

  double wz_sum = 1.0 / (static_cast<double>(beamspot.sigma_cm[2]) * beamspot.sigma_cm[2]);
  double bz = wz_sum * static_cast<double>(beamspot.xyz_cm[2]);

  const double w_rp = 1.0 / (static_cast<double>(track_sigma_rphi_cm) * track_sigma_rphi_cm);
  const double w_z = 1.0 / (static_cast<double>(track_sigma_z_cm) * track_sigma_z_cm);

  for (const auto& t : tracks) {
    // Perigee constraint: -x*sin(phi0) + y*cos(phi0) = eps
    const double s = std::sin(static_cast<double>(t.phi0));
    const double c = std::cos(static_cast<double>(t.phi0));
    const double nx = -s;
    const double ny = c;
    const double eps = static_cast<double>(t.eps);

    a11 += w_rp * nx * nx;
    a12 += w_rp * nx * ny;
    a22 += w_rp * ny * ny;
    bx += w_rp * nx * eps;
    by += w_rp * ny * eps;

    wz_sum += w_z;
    bz += w_z * static_cast<double>(t.z0);
  }

  const double det = a11 * a22 - a12 * a12;
  const double vx = (a22 * bx - a12 * by) / det;
  const double vy = (-a12 * bx + a11 * by) / det;
  const double vz = bz / wz_sum;

  double chi2 = 0.0;
  for (const auto& t : tracks) {
    const double s = std::sin(static_cast<double>(t.phi0));
    const double c = std::cos(static_cast<double>(t.phi0));
    const double r_rp = (-vx * s + vy * c - static_cast<double>(t.eps)) /
                        static_cast<double>(track_sigma_rphi_cm);
    const double r_z = (vz - static_cast<double>(t.z0)) /
                       static_cast<double>(track_sigma_z_cm);
    chi2 += r_rp * r_rp + r_z * r_z;
  }

  return VertexFitResult{
      .xyz_cm = {static_cast<float>(vx), static_cast<float>(vy), static_cast<float>(vz)},
      .sigma_cm = {static_cast<float>(std::sqrt(a22 / det)),
                   static_cast<float>(std::sqrt(a11 / det)),
                   static_cast<float>(std::sqrt(1.0 / wz_sum))},
      .chi2 = static_cast<float>(chi2),
      .ndf = static_cast<std::int32_t>(tracks.size() * 2),
      .num_tracks = static_cast<std::int32_t>(tracks.size()),
  };
}

}  // namespace delphipp::reco::tracking
