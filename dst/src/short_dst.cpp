// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/dst/short_dst.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace delphipp::dst {

namespace {

ParticleTag classify_charged(const TrackPidInput& in, float p_gev) noexcept {
  if (in.has_muon_hit && p_gev > 2.0f) {
    return ParticleTag::Muon;
  }
  if (p_gev > 1.0f && (in.ecal_energy_gev / p_gev) > 0.75f &&
      in.hcal_energy_gev < 1.0f) {
    return ParticleTag::Electron;
  }
  if (in.has_rich) {
    switch (in.rich_species) {
      case reco::rich::Species::Electron: return ParticleTag::Electron;
      case reco::rich::Species::Muon:     return ParticleTag::Muon;
      case reco::rich::Species::Pion:     return ParticleTag::Pion;
      case reco::rich::Species::Kaon:     return ParticleTag::Kaon;
      case reco::rich::Species::Proton:   return ParticleTag::Proton;
    }
  }
  return ParticleTag::Pion;
}

float tag_mass_gev(ParticleTag tag) noexcept {
  switch (tag) {
    case ParticleTag::Electron: return 0.000511f;
    case ParticleTag::Muon:     return 0.105658f;
    case ParticleTag::Pion:     return 0.139570f;
    case ParticleTag::Kaon:     return 0.493677f;
    case ParticleTag::Proton:   return 0.938272f;
    case ParticleTag::Photon:   return 0.0f;
    case ParticleTag::NeutralHadron: return 0.497611f;
    default:                    return 0.139570f;
  }
}

}  // namespace

DstEvent build_short_dst(std::int32_t run,
                         std::int32_t event,
                         float bz_gev_per_cm,
                         std::span<const TrackPidInput> charged_tracks,
                         std::span<const NeutralInput> neutrals,
                         const reco::tracking::BeamspotConstraint& beamspot,
                         bool is_bhabha) {
  DstEvent dst;
  dst.run = run;
  dst.event = event;
  dst.is_bhabha = is_bhabha;

  std::vector<PerigeeTrack> perigees;
  perigees.reserve(charged_tracks.size());
  for (const auto& ct : charged_tracks) {
    perigees.push_back(ct.track);
  }
  dst.primary_vertex = reco::tracking::fit_primary_vertex(perigees, beamspot);

  std::int32_t next_id = 1;
  float vis_e = 0.0f;
  float btag_sig = 0.0f;

  for (const auto& ct : charged_tracks) {
    const float p = ct.track.p_gev(bz_gev_per_cm);
    if (!std::isfinite(p) || p > 500.0f) {
      continue;
    }
    const float st = std::sin(ct.track.theta);
    const float ct_th = std::cos(ct.track.theta);
    const float px = p * st * std::cos(ct.track.phi0);
    const float py = p * st * std::sin(ct.track.phi0);
    const float pz = p * ct_th;
    const ParticleTag tag = classify_charged(ct, p);
    const float m = tag_mass_gev(tag);
    const float energy = std::hypot(p, m);

    // Impact parameter w.r.t. primary vertex
    const float ip_rphi = ct.track.eps -
                          (-dst.primary_vertex.xyz_cm[0] * std::sin(ct.track.phi0) +
                           dst.primary_vertex.xyz_cm[1] * std::cos(ct.track.phi0));
    const float ip_z = ct.track.z0 - dst.primary_vertex.xyz_cm[2];
    btag_sig += std::fabs(ip_rphi) / 0.005f;

    dst.particles.push_back(DstParticle{
        .id = next_id++,
        .charge = ct.charge,
        .px = px,
        .py = py,
        .pz = pz,
        .energy = energy,
        .impact_rphi_cm = ip_rphi,
        .impact_z_cm = ip_z,
        .ecal_energy_gev = ct.ecal_energy_gev,
        .hcal_energy_gev = ct.hcal_energy_gev,
        .tpc_dedx = ct.tpc_dedx,
        .muon_tag = (tag == ParticleTag::Muon),
        .pid = tag,
    });
    vis_e += energy;
  }

  for (const auto& nu : neutrals) {
    const float r = std::sqrt(nu.xyz_cm[0] * nu.xyz_cm[0] +
                              nu.xyz_cm[1] * nu.xyz_cm[1] +
                              nu.xyz_cm[2] * nu.xyz_cm[2]);
    if (r <= 1e-6f || nu.energy_gev <= 0.0f) continue;
    const ParticleTag tag = nu.is_ecal ? ParticleTag::Photon : ParticleTag::NeutralHadron;
    dst.particles.push_back(DstParticle{
        .id = next_id++,
        .charge = 0,
        .px = nu.energy_gev * nu.xyz_cm[0] / r,
        .py = nu.energy_gev * nu.xyz_cm[1] / r,
        .pz = nu.energy_gev * nu.xyz_cm[2] / r,
        .energy = nu.energy_gev,
        .ecal_energy_gev = nu.is_ecal ? nu.energy_gev : 0.0f,
        .hcal_energy_gev = nu.is_ecal ? 0.0f : nu.energy_gev,
        .pid = tag,
    });
    vis_e += nu.energy_gev;
  }

  dst.visible_energy_gev = vis_e;
  dst.btag_significance_sum = btag_sig;
  return dst;
}

std::string DstEvent::to_json() const {
  std::ostringstream os;
  os << std::setprecision(6);
  os << "{\"run\":" << run << ",\"event\":" << event
     << ",\"visible_energy_gev\":" << visible_energy_gev
     << ",\"is_bhabha\":" << (is_bhabha ? "true" : "false")
     << ",\"primary_vertex\":[" << primary_vertex.xyz_cm[0] << ','
     << primary_vertex.xyz_cm[1] << ',' << primary_vertex.xyz_cm[2] << ']'
     << ",\"particles\":[";
  for (std::size_t i = 0; i < particles.size(); ++i) {
    if (i > 0) os << ',';
    const auto& p = particles[i];
    os << "{\"id\":" << p.id << ",\"charge\":" << p.charge
       << ",\"pid\":" << static_cast<std::int32_t>(p.pid)
       << ",\"px\":" << p.px << ",\"py\":" << p.py << ",\"pz\":" << p.pz
       << ",\"E\":" << p.energy << '}';
  }
  os << "]}";
  return os.str();
}

}  // namespace delphipp::dst
