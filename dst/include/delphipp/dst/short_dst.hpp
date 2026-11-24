// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 DELPHI Short/Long DST data model and builder
// (simana_v94d/src/car/pxdst33.car, phdstxx.car, shortdst.car, muflag.car,
// elepha32.car, aabtagxx.car).

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "delphipp/geometry.hpp"
#include "delphipp/helix.hpp"
#include "delphipp/reco/rich.hpp"
#include "delphipp/reco/vertex_fitter.hpp"
#include "delphipp/tanagra_store.hpp"

namespace delphipp::dst {

enum class ParticleTag : std::int32_t {
  Unknown = 0,
  Electron = 11,
  Muon = 13,
  Pion = 211,
  Kaon = 321,
  Proton = 2212,
  Photon = 22,
  NeutralHadron = 130,
};

struct DstParticle {
  std::int32_t id = 0;
  std::int32_t charge = 0;  // -1, 0, +1
  float px = 0.0f;
  float py = 0.0f;
  float pz = 0.0f;
  float energy = 0.0f;
  float impact_rphi_cm = 0.0f;
  float impact_z_cm = 0.0f;
  float ecal_energy_gev = 0.0f;
  float hcal_energy_gev = 0.0f;
  float tpc_dedx = 0.0f;
  bool muon_tag = false;
  ParticleTag pid = ParticleTag::Unknown;
};

struct DstEvent {
  std::int32_t run = 0;
  std::int32_t event = 0;
  float ecm_gev = 91.2f;
  reco::tracking::VertexFitResult primary_vertex{};
  std::vector<DstParticle> particles;
  float visible_energy_gev = 0.0f;
  float btag_significance_sum = 0.0f;
  bool is_bhabha = false;

  [[nodiscard]] std::string to_json() const;
};

struct TrackPidInput {
  PerigeeTrack track;
  std::int32_t charge = +1;
  float ecal_energy_gev = 0.0f;
  float hcal_energy_gev = 0.0f;
  float tpc_dedx = 1.0f;
  bool has_muon_hit = false;
  reco::rich::Species rich_species = reco::rich::Species::Pion;
  bool has_rich = false;
};

struct NeutralInput {
  Vec3 xyz_cm{};
  float energy_gev = 0.0f;
  bool is_ecal = true;  // true -> photon, false -> neutral hadron
};

// Build a complete Short-DST event from reconstructed tracks, neutrals, field, and beamspot.
[[nodiscard]] DstEvent build_short_dst(
    std::int32_t run,
    std::int32_t event,
    float bz_gev_per_cm,
    std::span<const TrackPidInput> charged_tracks,
    std::span<const NeutralInput> neutrals,
    const reco::tracking::BeamspotConstraint& beamspot = {},
    bool is_bhabha = false);

}  // namespace delphipp::dst
