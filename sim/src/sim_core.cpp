// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/sim/sim_core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace delphipp::sim {

namespace {

constexpr double kPi = std::numbers::pi_v<double>;
constexpr double kTwoPi = 2.0 * kPi;

double normalize_phi(double phi) noexcept {
  while (phi < 0.0) {
    phi += kTwoPi;
  }
  while (phi >= kTwoPi) {
    phi -= kTwoPi;
  }
  return phi;
}

double resolve_mass_gev(std::int32_t pdg, float explicit_mass) noexcept {
  if (std::isfinite(explicit_mass) && explicit_mass > 0.0f) {
    return static_cast<double>(explicit_mass);
  }
  const std::int64_t pdg64 = static_cast<std::int64_t>(pdg);
  const std::int64_t apdg = (pdg64 < 0) ? -pdg64 : pdg64;
  switch (apdg) {
    case 11:   return 0.00051099895;
    case 13:   return 0.1056583755;
    case 211:  return 0.13957039;
    case 321:  return 0.493677;
    case 2212: return 0.938272088;
    default:   return 0.13957039;
  }
}

struct Vec3d {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3d normalize_vec(const Vec3d& v) noexcept {
  const double n = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  if (n <= 1e-12) {
    return {1.0, 0.0, 0.0};
  }
  return {v.x / n, v.y / n, v.z / n};
}

// Lorentz force derivative d(dir)/ds = (q / p) * (B(pos) x dir),
// matching delphipp::helix curvature sign k = q * Bz / pt.
Vec3d lorentz_dudt(const Vec3d& pos,
                   const Vec3d& dir,
                   double q_over_p,
                   const MagneticField& bfield) noexcept {
  const FieldVector bv = bfield.evaluate(
      {static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z)});
  const double bx = static_cast<double>(bv.bx);
  const double by = static_cast<double>(bv.by);
  const double bz = static_cast<double>(bv.bz);
  return {
      q_over_p * (by * dir.z - bz * dir.y),
      q_over_p * (bz * dir.x - bx * dir.z),
      q_over_p * (bx * dir.y - by * dir.x),
  };
}

void rk4_step(Vec3d& pos,
              Vec3d& dir,
              double ds,
              double q_over_p,
              const MagneticField& bfield) noexcept {
  const Vec3d k1_r = dir;
  const Vec3d k1_u = lorentz_dudt(pos, dir, q_over_p, bfield);

  const double h2 = 0.5 * ds;
  const Vec3d p2{pos.x + h2 * k1_r.x, pos.y + h2 * k1_r.y, pos.z + h2 * k1_r.z};
  const Vec3d u2 = normalize_vec({dir.x + h2 * k1_u.x, dir.y + h2 * k1_u.y, dir.z + h2 * k1_u.z});
  const Vec3d k2_r = u2;
  const Vec3d k2_u = lorentz_dudt(p2, u2, q_over_p, bfield);

  const Vec3d p3{pos.x + h2 * k2_r.x, pos.y + h2 * k2_r.y, pos.z + h2 * k2_r.z};
  const Vec3d u3 = normalize_vec({dir.x + h2 * k2_u.x, dir.y + h2 * k2_u.y, dir.z + h2 * k2_u.z});
  const Vec3d k3_r = u3;
  const Vec3d k3_u = lorentz_dudt(p3, u3, q_over_p, bfield);

  const Vec3d p4{pos.x + ds * k3_r.x, pos.y + ds * k3_r.y, pos.z + ds * k3_r.z};
  const Vec3d u4 = normalize_vec({dir.x + ds * k3_u.x, dir.y + ds * k3_u.y, dir.z + ds * k3_u.z});
  const Vec3d k4_r = u4;
  const Vec3d k4_u = lorentz_dudt(p4, u4, q_over_p, bfield);

  const double c6 = ds / 6.0;
  pos.x += c6 * (k1_r.x + 2.0 * k2_r.x + 2.0 * k3_r.x + k4_r.x);
  pos.y += c6 * (k1_r.y + 2.0 * k2_r.y + 2.0 * k3_r.y + k4_r.y);
  pos.z += c6 * (k1_r.z + 2.0 * k2_r.z + 2.0 * k3_r.z + k4_r.z);

  dir = normalize_vec({
      dir.x + c6 * (k1_u.x + 2.0 * k2_u.x + 2.0 * k3_u.x + k4_u.x),
      dir.y + c6 * (k1_u.y + 2.0 * k2_u.y + 2.0 * k3_u.y + k4_u.y),
      dir.z + c6 * (k1_u.z + 2.0 * k2_u.z + 2.0 * k3_u.z + k4_u.z),
  });
}

double highland_rms_rad(double p_gev, double beta, double x_over_x0) noexcept {
  if (p_gev <= 1e-4 || x_over_x0 <= 1e-9) {
    return 0.0;
  }
  const double b = std::clamp(beta, 0.05, 1.0);
  const double sq = std::sqrt(x_over_x0);
  const double corr = std::max(0.1, 1.0 + 0.038 * std::log(x_over_x0));
  return (0.0136 / (b * p_gev)) * sq * corr;
}

double bethe_bloch_dedx_mev(double p_gev,
                            double e_gev,
                            double mass_gev,
                            double abs_charge,
                            bool is_silicon,
                            double eff_thickness_cm) noexcept {
  const double beta = std::clamp(p_gev / std::max(e_gev, 1e-6), 0.05, 0.9999999);
  const double m_mev = std::max(mass_gev * 1000.0, 0.51099895);
  const double p_mev = p_gev * 1000.0;
  const double e_mev = e_gev * 1000.0;
  const double bg = std::max(p_mev / m_mev, 0.05);
  const double gamma = std::max(e_mev / m_mev, 1.0);

  constexpr double kBethe = 0.307075;  // MeV * cm^2 / g
  constexpr double kMeMev = 0.51099895;
  const double z_over_a = is_silicon ? 0.4985 : 0.450;
  const double i_mev = is_silicon ? 173.0e-6 : 188.0e-6;
  const double rho = is_silicon ? 2.329 : 0.0017;

  const double me_ratio = kMeMev / m_mev;
  const double w_max = (2.0 * kMeMev * bg * bg) /
                       (1.0 + 2.0 * gamma * me_ratio + me_ratio * me_ratio);
  const double arg = std::max(1.1, (2.0 * kMeMev * bg * bg * w_max) / (i_mev * i_mev));
  const double stopping =
      kBethe * abs_charge * abs_charge * z_over_a / (beta * beta) *
      (0.5 * std::log(arg) - beta * beta);
  const double clamped_stopping = std::clamp(stopping, 1.0, 50.0);
  return clamped_stopping * rho * std::max(eff_thickness_cm, 0.001);
}

void apply_highland_ms(Vec3d& dir,
                       double theta0,
                       std::mt19937_64& rng,
                       std::normal_distribution<double>& gauss) noexcept {
  if (theta0 <= 1e-8) {
    return;
  }
  const Vec3d ref = (std::fabs(dir.z) < 0.9) ? Vec3d{0.0, 0.0, 1.0} : Vec3d{1.0, 0.0, 0.0};
  const Vec3d e1 = normalize_vec({
      dir.y * ref.z - dir.z * ref.y,
      dir.z * ref.x - dir.x * ref.z,
      dir.x * ref.y - dir.y * ref.x,
  });
  const Vec3d e2 = normalize_vec({
      dir.y * e1.z - dir.z * e1.y,
      dir.z * e1.x - dir.x * e1.z,
      dir.x * e1.y - dir.y * e1.x,
  });
  const double th1 = theta0 * gauss(rng);
  const double th2 = theta0 * gauss(rng);
  dir = normalize_vec({
      dir.x + th1 * e1.x + th2 * e2.x,
      dir.y + th1 * e1.y + th2 * e2.y,
      dir.z + th1 * e1.z + th2 * e2.z,
  });
}

struct BarrelSurface {
  const char* name;
  double r_cm;
  double max_abs_z_cm;
  double x0_normal;
  double thickness_cm;
  bool is_silicon;
};

constexpr std::array<BarrelSurface, 26> kBarrelSurfaces = {{
    {"VD",  6.30,  18.0, 0.0032, 0.030, true},
    {"VD",  9.00,  18.0, 0.0032, 0.030, true},
    {"VD",  10.90, 18.0, 0.0032, 0.030, true},
    {"ID",  12.26, 50.0, 0.0015, 1.80,  false},
    {"ID",  14.08, 50.0, 0.0015, 1.80,  false},
    {"ID",  15.91, 50.0, 0.0015, 1.80,  false},
    {"ID",  17.73, 50.0, 0.0015, 1.80,  false},
    {"ID",  19.56, 50.0, 0.0015, 1.80,  false},
    {"ID",  21.39, 50.0, 0.0015, 1.80,  false},
    {"STW", 23.20, 50.0, 0.0016, 1.20,  false},
    {"STW", 24.40, 50.0, 0.0016, 1.20,  false},
    {"STW", 25.60, 50.0, 0.0016, 1.20,  false},
    {"STW", 26.80, 50.0, 0.0016, 1.20,  false},
    {"STW", 28.00, 50.0, 0.0016, 1.20,  false},
    {"TPC", 35.00, 134.0, 0.0005, 5.00, false},
    {"TPC", 45.00, 134.0, 0.0005, 5.00, false},
    {"TPC", 55.00, 134.0, 0.0005, 5.00, false},
    {"TPC", 65.00, 134.0, 0.0005, 5.00, false},
    {"TPC", 80.00, 134.0, 0.0005, 5.00, false},
    {"TPC", 95.00, 134.0, 0.0005, 5.00, false},
    {"TPC", 110.00, 134.0, 0.0005, 5.00, false},
    {"OD",  197.00, 235.0, 0.0030, 2.25, false},
    {"OD",  199.25, 235.0, 0.0030, 2.25, false},
    {"OD",  201.50, 235.0, 0.0030, 2.25, false},
    {"OD",  203.75, 235.0, 0.0030, 2.25, false},
    {"OD",  206.00, 235.0, 0.0030, 2.25, false},
}};

struct ForwardSurface {
  const char* name;
  double z_cm;
  double r_min_cm;
  double r_max_cm;
  double x0_normal;
  double thickness_cm;
};

constexpr std::array<ForwardSurface, 4> kForwardSurfaces = {{
    {"FCA",  160.0, 25.0, 155.0, 0.015, 4.0},
    {"FCA", -160.0, 25.0, 155.0, 0.015, 4.0},
    {"FCB",  275.0, 45.0, 200.0, 0.020, 5.0},
    {"FCB", -275.0, 45.0, 200.0, 0.020, 5.0},
}};

}  // namespace

PrimaryVertexSmearer::PrimaryVertexSmearer(std::uint64_t seed,
                                           BeamspotProfile profile) noexcept
    : profile_(profile), rng_(seed) {}

PrimaryVertexSmearer::PrimaryVertexSmearer(std::uint64_t seed,
                                           float sigma_x_cm,
                                           float sigma_y_cm,
                                           float sigma_z_cm,
                                           Vec3 nominal_cm) noexcept
    : profile_{
          .nominal_cm = nominal_cm,
          .sigma_x_cm = sigma_x_cm,
          .sigma_y_cm = sigma_y_cm,
          .sigma_z_cm = sigma_z_cm,
      },
      rng_(seed) {}

Vec3 PrimaryVertexSmearer::sample() noexcept {
  std::normal_distribution<double> gauss(0.0, 1.0);
  auto clean_sigma = [](float s) noexcept {
    return (std::isfinite(s) && s > 0.0f) ? static_cast<double>(s) : 0.0;
  };
  auto clean_nom = [](float c) noexcept {
    return std::isfinite(c) ? static_cast<double>(c) : 0.0;
  };
  const double vx = clean_nom(profile_.nominal_cm[0]) + clean_sigma(profile_.sigma_x_cm) * gauss(rng_);
  const double vy = clean_nom(profile_.nominal_cm[1]) + clean_sigma(profile_.sigma_y_cm) * gauss(rng_);
  const double vz = clean_nom(profile_.nominal_cm[2]) + clean_sigma(profile_.sigma_z_cm) * gauss(rng_);
  return {static_cast<float>(vx), static_cast<float>(vy), static_cast<float>(vz)};
}

std::vector<SteppedHit> step_particle_in_bfield(const SimParticle& particle,
                                                 const MagneticField& bfield,
                                                 std::uint64_t seed) {
  std::vector<SteppedHit> hits;
  if (particle.charge == 0) {
    return hits;
  }
  if (!std::isfinite(particle.px) || !std::isfinite(particle.py) ||
      !std::isfinite(particle.pz) || !std::isfinite(particle.vertex_cm[0]) ||
      !std::isfinite(particle.vertex_cm[1]) || !std::isfinite(particle.vertex_cm[2])) {
    return hits;
  }

  double p_gev = std::sqrt(static_cast<double>(particle.px) * particle.px +
                           static_cast<double>(particle.py) * particle.py +
                           static_cast<double>(particle.pz) * particle.pz);
  if (!std::isfinite(p_gev) || p_gev < 1e-4) {
    return hits;
  }

  const double mass_gev = resolve_mass_gev(particle.pdg, particle.mass);
  double e_gev = std::sqrt(p_gev * p_gev + mass_gev * mass_gev);
  const double q = static_cast<double>(particle.charge);
  const double abs_q = std::fabs(q);

  Vec3d pos{
      static_cast<double>(particle.vertex_cm[0]),
      static_cast<double>(particle.vertex_cm[1]),
      static_cast<double>(particle.vertex_cm[2]),
  };
  Vec3d dir = normalize_vec({
      static_cast<double>(particle.px),
      static_cast<double>(particle.py),
      static_cast<double>(particle.pz),
  });

  std::mt19937_64 rng(seed);
  std::normal_distribution<double> gauss(0.0, 1.0);

  std::array<bool, kBarrelSurfaces.size()> barrel_visited{};
  std::array<bool, kForwardSurfaces.size()> fwd_visited{};

  double s_total = 0.0;
  constexpr double kMaxPathCm = 600.0;
  constexpr std::size_t kMaxSteps = 1200;

  for (std::size_t step = 0; step < kMaxSteps && s_total < kMaxPathCm && p_gev >= 0.02; ++step) {
    const double r_prev = std::hypot(pos.x, pos.y);
    if (r_prev > 215.0 || std::fabs(pos.z) > 290.0) {
      break;
    }

    const double ds_nom = (r_prev < 32.0) ? 0.5 : 1.0;
    const double q_over_p = q / p_gev;

    Vec3d pos_next = pos;
    Vec3d dir_next = dir;
    rk4_step(pos_next, dir_next, ds_nom, q_over_p, bfield);
    const double r_next = std::hypot(pos_next.x, pos_next.y);

    double best_t = 2.0;
    bool is_barrel_hit = false;
    std::size_t hit_idx = 0;

    if (r_next > r_prev) {
      for (std::size_t i = 0; i < kBarrelSurfaces.size(); ++i) {
        if (barrel_visited[i]) {
          continue;
        }
        const double rt = kBarrelSurfaces[i].r_cm;
        if (r_prev < rt && r_next >= rt) {
          const double t = std::clamp((rt - r_prev) / (r_next - r_prev), 1e-4, 1.0);
          if (t < best_t) {
            best_t = t;
            is_barrel_hit = true;
            hit_idx = i;
          }
        }
      }
    }

    for (std::size_t i = 0; i < kForwardSurfaces.size(); ++i) {
      if (fwd_visited[i]) {
        continue;
      }
      const double zt = kForwardSurfaces[i].z_cm;
      const bool crossed_pos = (zt > 0.0 && pos.z < zt && pos_next.z >= zt);
      const bool crossed_neg = (zt < 0.0 && pos.z > zt && pos_next.z <= zt);
      if (crossed_pos || crossed_neg) {
        const double dz = pos_next.z - pos.z;
        if (std::fabs(dz) > 1e-12) {
          const double t = std::clamp((zt - pos.z) / dz, 1e-4, 1.0);
          if (t < best_t) {
            best_t = t;
            is_barrel_hit = false;
            hit_idx = i;
          }
        }
      }
    }

    if (best_t <= 1.0) {
      const double ds_sub = best_t * ds_nom;
      rk4_step(pos, dir, ds_sub, q_over_p, bfield);
      s_total += std::max(ds_sub, 1e-4);

      if (is_barrel_hit) {
        barrel_visited[hit_idx] = true;
        const auto& surf = kBarrelSurfaces[hit_idx];
        if (std::fabs(pos.z) <= surf.max_abs_z_cm) {
          const double sin_th = std::max(0.15, std::hypot(dir.x, dir.y));
          const double path_x0 = surf.x0_normal / sin_th;
          const double eff_thick = surf.thickness_cm / sin_th;
          const double dedx_mev = bethe_bloch_dedx_mev(
              p_gev, e_gev, mass_gev, abs_q, surf.is_silicon, eff_thick);

          hits.push_back(SteppedHit{
              .det_name = surf.name,
              .r_cm = static_cast<float>(surf.r_cm),
              .phi_rad = static_cast<float>(normalize_phi(std::atan2(pos.y, pos.x))),
              .z_cm = static_cast<float>(pos.z),
              .path_x0 = static_cast<float>(path_x0),
              .dedx_mev = static_cast<float>(dedx_mev),
          });

          const double beta = std::clamp(p_gev / std::max(e_gev, 1e-6), 0.05, 1.0);
          const double th0 = highland_rms_rad(p_gev, beta, path_x0);
          apply_highland_ms(dir, th0, rng, gauss);

          e_gev = std::max(mass_gev + 0.001, e_gev - dedx_mev * 1e-3);
          p_gev = std::sqrt(std::max(0.0, e_gev * e_gev - mass_gev * mass_gev));
        }
      } else {
        fwd_visited[hit_idx] = true;
        const auto& surf = kForwardSurfaces[hit_idx];
        const double r_hit = std::hypot(pos.x, pos.y);
        if (r_hit >= surf.r_min_cm && r_hit <= surf.r_max_cm) {
          const double cos_th = std::max(0.15, std::fabs(dir.z));
          const double path_x0 = surf.x0_normal / cos_th;
          const double eff_thick = surf.thickness_cm / cos_th;
          const double dedx_mev = bethe_bloch_dedx_mev(
              p_gev, e_gev, mass_gev, abs_q, false, eff_thick);

          hits.push_back(SteppedHit{
              .det_name = surf.name,
              .r_cm = static_cast<float>(r_hit),
              .phi_rad = static_cast<float>(normalize_phi(std::atan2(pos.y, pos.x))),
              .z_cm = static_cast<float>(surf.z_cm),
              .path_x0 = static_cast<float>(path_x0),
              .dedx_mev = static_cast<float>(dedx_mev),
          });

          const double beta = std::clamp(p_gev / std::max(e_gev, 1e-6), 0.05, 1.0);
          const double th0 = highland_rms_rad(p_gev, beta, path_x0);
          apply_highland_ms(dir, th0, rng, gauss);

          e_gev = std::max(mass_gev + 0.001, e_gev - dedx_mev * 1e-3);
          p_gev = std::sqrt(std::max(0.0, e_gev * e_gev - mass_gev * mass_gev));
        }
      }
    } else {
      pos = pos_next;
      dir = dir_next;
      s_total += ds_nom;
    }
  }

  return hits;
}

}  // namespace delphipp::sim
