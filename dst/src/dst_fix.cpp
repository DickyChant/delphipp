// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/dst/dst_fix.hpp"

#include <cmath>
#include <vector>

namespace delphipp::dst {

void apply_dst_fixes(DstEvent& ev, const DstFixConfig& cfg) {
  std::vector<DstParticle> kept;
  kept.reserve(ev.particles.size());
  float vis_e = 0.0f;

  for (auto p : ev.particles) {
    if (p.charge == 0) {
      const float scale = (p.pid == ParticleTag::NeutralHadron) ? cfg.hcal_scale
                                                                : cfg.ecal_scale;
      p.px *= scale;
      p.py *= scale;
      p.pz *= scale;
      p.energy *= scale;
      p.ecal_energy_gev *= cfg.ecal_scale;
      p.hcal_energy_gev *= cfg.hcal_scale;
      if (!std::isfinite(p.energy) || p.energy < cfg.min_neutral_energy_gev) {
        continue;
      }
    } else {
      const float sign_z = (p.pz >= 0.0f) ? 1.0f : -1.0f;
      p.impact_z_cm += sign_z * cfg.tpc_z_shift_cm;
      p.ecal_energy_gev *= cfg.ecal_scale;
      p.hcal_energy_gev *= cfg.hcal_scale;
    }
    if (std::isfinite(p.energy) && p.energy > 0.0f) {
      vis_e += p.energy;
    }
    kept.push_back(p);
  }

  ev.particles = std::move(kept);
  ev.visible_energy_gev = vis_e;
}

}  // namespace delphipp::dst
