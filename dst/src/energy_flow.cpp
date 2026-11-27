// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/dst/energy_flow.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace delphipp::dst {

namespace {

double durham_yij(const Jet& a, const Jet& b, double evis2) noexcept {
  if (evis2 <= 1e-12) {
    return 0.0;
  }
  const double pa = std::sqrt(static_cast<double>(a.px) * a.px +
                              static_cast<double>(a.py) * a.py +
                              static_cast<double>(a.pz) * a.pz);
  const double pb = std::sqrt(static_cast<double>(b.px) * b.px +
                              static_cast<double>(b.py) * b.py +
                              static_cast<double>(b.pz) * b.pz);
  if (pa <= 1e-9 || pb <= 1e-9) {
    return 0.0;
  }
  const double dot = static_cast<double>(a.px) * b.px +
                     static_cast<double>(a.py) * b.py +
                     static_cast<double>(a.pz) * b.pz;
  const double cos_th = std::clamp(dot / (pa * pb), -1.0, 1.0);
  const double min_e2 = std::min(static_cast<double>(a.energy) * a.energy,
                                 static_cast<double>(b.energy) * b.energy);
  return 2.0 * min_e2 * (1.0 - cos_th) / evis2;
}

std::vector<Jet> init_proto_jets(std::span<const DstParticle> particles,
                                 double& total_e) {
  std::vector<Jet> jets;
  jets.reserve(particles.size());
  total_e = 0.0;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const auto& p = particles[i];
    if (!std::isfinite(p.energy) || p.energy <= 0.0f) {
      continue;
    }
    total_e += static_cast<double>(p.energy);
    jets.push_back(Jet{
        .px = p.px,
        .py = p.py,
        .pz = p.pz,
        .energy = p.energy,
        .constituent_indices = {i},
    });
  }
  return jets;
}

void merge_pair(std::vector<Jet>& jets, std::size_t i, std::size_t j) {
  jets[i].px += jets[j].px;
  jets[i].py += jets[j].py;
  jets[i].pz += jets[j].pz;
  jets[i].energy += jets[j].energy;
  jets[i].constituent_indices.insert(
      jets[i].constituent_indices.end(),
      jets[j].constituent_indices.begin(),
      jets[j].constituent_indices.end());
  jets.erase(jets.begin() + static_cast<std::ptrdiff_t>(j));
}

}  // namespace

std::vector<Jet> cluster_durham_jets(std::span<const DstParticle> particles,
                                     float y_cut,
                                     std::size_t min_jets) {
  double total_e = 0.0;
  auto jets = init_proto_jets(particles, total_e);
  if (jets.size() <= min_jets || total_e <= 1e-6) {
    return jets;
  }
  const double evis2 = total_e * total_e;
  const double cut = static_cast<double>(std::max(0.0f, y_cut));

  while (jets.size() > std::max<std::size_t>(1, min_jets)) {
    double best_y = std::numeric_limits<double>::infinity();
    std::size_t bi = 0, bj = 1;
    for (std::size_t i = 0; i < jets.size(); ++i) {
      for (std::size_t j = i + 1; j < jets.size(); ++j) {
        const double y = durham_yij(jets[i], jets[j], evis2);
        if (y < best_y) {
          best_y = y;
          bi = i;
          bj = j;
        }
      }
    }
    if (best_y > cut) {
      break;
    }
    merge_pair(jets, bi, bj);
  }

  std::sort(jets.begin(), jets.end(), [](const Jet& a, const Jet& b) {
    return a.energy > b.energy;
  });
  return jets;
}

EventShape compute_event_shape(std::span<const DstParticle> particles) {
  EventShape shape{};
  if (particles.empty()) {
    return shape;
  }

  double sum_e = 0.0, sum_px = 0.0, sum_py = 0.0, sum_pz = 0.0, sum_pmag = 0.0;
  struct Mom3 {
    double x, y, z;
  };
  std::vector<Mom3> p3;
  p3.reserve(particles.size());

  for (const auto& pt : particles) {
    if (!std::isfinite(pt.px) || !std::isfinite(pt.py) ||
        !std::isfinite(pt.pz) || !std::isfinite(pt.energy) ||
        pt.energy <= 0.0f) {
      continue;
    }
    const double px = static_cast<double>(pt.px);
    const double py = static_cast<double>(pt.py);
    const double pz = static_cast<double>(pt.pz);
    sum_e += static_cast<double>(pt.energy);
    sum_px += px;
    sum_py += py;
    sum_pz += pz;
    sum_pmag += std::sqrt(px * px + py * py + pz * pz);
    p3.push_back({px, py, pz});
  }

  if (p3.empty() || sum_pmag <= 1e-9) {
    return shape;
  }

  const double m2 = sum_e * sum_e - (sum_px * sum_px + sum_py * sum_py + sum_pz * sum_pz);
  shape.visible_mass_gev = static_cast<float>(std::sqrt(std::max(0.0, m2)));
  shape.missing_pt_gev = static_cast<float>(std::hypot(sum_px, sum_py));

  // Evaluate candidate thrust axes from particle directions & all pairwise sums
  double best_t = 0.0;
  Mom3 best_n{0.0, 0.0, 1.0};
  auto eval_axis = [&](double nx, double ny, double nz) {
    const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (norm <= 1e-12) {
      return;
    }
    nx /= norm;
    ny /= norm;
    nz /= norm;
    // One Brandt-Dahmen fixed-point iteration from this seed
    for (int iter = 0; iter < 3; ++iter) {
      double sx = 0.0, sy = 0.0, sz = 0.0;
      for (const auto& v : p3) {
        const double sgn = (v.x * nx + v.y * ny + v.z * nz >= 0.0) ? 1.0 : -1.0;
        sx += sgn * v.x;
        sy += sgn * v.y;
        sz += sgn * v.z;
      }
      const double sn = std::sqrt(sx * sx + sy * sy + sz * sz);
      if (sn <= 1e-12) {
        break;
      }
      nx = sx / sn;
      ny = sy / sn;
      nz = sz / sn;
    }
    double num = 0.0;
    for (const auto& v : p3) {
      num += std::fabs(v.x * nx + v.y * ny + v.z * nz);
    }
    const double t = num / sum_pmag;
    if (t > best_t) {
      best_t = t;
      best_n = {nx, ny, nz};
    }
  };

  for (std::size_t i = 0; i < p3.size(); ++i) {
    eval_axis(p3[i].x, p3[i].y, p3[i].z);
    for (std::size_t j = i + 1; j < std::min<std::size_t>(p3.size(), 24); ++j) {
      eval_axis(p3[i].x + p3[j].x, p3[i].y + p3[j].y, p3[i].z + p3[j].z);
    }
  }
  shape.thrust = static_cast<float>(std::clamp(best_t, 0.5, 1.0));
  shape.thrust_axis = {
      static_cast<float>(best_n.x),
      static_cast<float>(best_n.y),
      static_cast<float>(best_n.z),
  };

  // Compute y_32 (Durham y_ij at transition from 3 jets -> 2 jets)
  double total_e = 0.0;
  auto jets = init_proto_jets(particles, total_e);
  if (jets.size() >= 3 && total_e > 1e-6) {
    const double evis2 = total_e * total_e;
    while (jets.size() > 2) {
      double best_y = std::numeric_limits<double>::infinity();
      std::size_t bi = 0, bj = 1;
      for (std::size_t i = 0; i < jets.size(); ++i) {
        for (std::size_t j = i + 1; j < jets.size(); ++j) {
          const double y = durham_yij(jets[i], jets[j], evis2);
          if (y < best_y) {
            best_y = y;
            bi = i;
            bj = j;
          }
        }
      }
      if (jets.size() == 3) {
        shape.y_cut_3to2 = static_cast<float>(best_y);
      }
      merge_pair(jets, bi, bj);
    }
  }

  return shape;
}

}  // namespace delphipp::dst
