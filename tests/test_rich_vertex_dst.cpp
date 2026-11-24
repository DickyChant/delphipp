// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

#include "delphipp/dst/short_dst.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/reco/rich.hpp"
#include "delphipp/reco/vertex_fitter.hpp"
#include "delphipp/tanagra_store.hpp"

int main() {
  using namespace delphipp;

  // 1. RICH Kaon identification in liquid C6F14 at p = 2.0 GeV/c
  const auto th_kaon = reco::rich::expected_cherenkov_angle(
      2.0f, reco::rich::kSpeciesMassGeV[static_cast<std::size_t>(reco::rich::Species::Kaon)],
      reco::rich::kRefractiveIndexLiquid);
  assert(th_kaon.has_value());
  reco::rich::RichRingMeasurement ring{
      .momentum_gev = 2.0f,
      .refractive_index = reco::rich::kRefractiveIndexLiquid,
      .measured_theta_rad = *th_kaon,
      .sigma_theta_rad = 0.005f,
      .num_photons = 12,
  };
  const auto rich_res = reco::rich::evaluate_hypotheses(ring);
  assert(rich_res.best_species == reco::rich::Species::Kaon);

  TanagraStore store;
  assert(reco::rich::write_ter_bank(store, ring, rich_res) == 1);

  // 2. Primary vertex fitter recovering shifted vertex (0.02, -0.01, 0.35) cm
  const float vx_true = 0.02f, vy_true = -0.01f, vz_true = 0.35f;
  std::vector<PerigeeTrack> vtx_tracks;
  for (int i = 0; i < 6; ++i) {
    const float phi = static_cast<float>(i) * (std::numbers::pi_v<float> / 3.0f);
    const float eps = -vx_true * std::sin(phi) + vy_true * std::cos(phi);
    vtx_tracks.push_back(PerigeeTrack{.eps = eps, .z0 = vz_true, .theta = 1.4f, .phi0 = phi, .inv_r = 0.0005f});
  }
  const auto vtx = reco::tracking::fit_primary_vertex(
      vtx_tracks, reco::tracking::BeamspotConstraint{.xyz_cm = {0.0f, 0.0f, 0.0f}, .sigma_cm = {0.1f, 0.1f, 1.0f}});
  assert(std::fabs(vtx.xyz_cm[0] - vx_true) < 0.005f);
  assert(std::fabs(vtx.xyz_cm[1] - vy_true) < 0.005f);
  assert(std::fabs(vtx.xyz_cm[2] - vz_true) < 0.01f);

  // 3. Full Short-DST builder on Z -> mu+ mu- (45.6 GeV each)
  MagneticField bfield;
  const float bz = bfield.evaluate({0.0f, 0.0f, 0.0f}).bz;
  const float inv_r_45 = bz / 45.6f;

  std::vector<dst::TrackPidInput> tracks = {
      {.track = {.eps = 0.0f, .z0 = 0.0f, .theta = std::numbers::pi_v<float> * 0.5f, .phi0 = 0.2f, .inv_r = inv_r_45},
       .charge = +1, .ecal_energy_gev = 0.2f, .hcal_energy_gev = 1.5f, .has_muon_hit = true},
      {.track = {.eps = 0.0f, .z0 = 0.0f, .theta = std::numbers::pi_v<float> * 0.5f, .phi0 = std::numbers::pi_v<float> + 0.2f, .inv_r = -inv_r_45},
       .charge = -1, .ecal_energy_gev = 0.2f, .hcal_energy_gev = 1.5f, .has_muon_hit = true},
  };
  std::vector<dst::NeutralInput> neutrals = {
      {.xyz_cm = {200.0f, 0.0f, 50.0f}, .energy_gev = 0.5f, .is_ecal = true},
  };

  const auto dst_ev = dst::build_short_dst(47014, 1, bz, tracks, neutrals);
  assert(dst_ev.particles.size() == 3);
  assert(dst_ev.particles[0].pid == dst::ParticleTag::Muon);
  assert(dst_ev.particles[1].pid == dst::ParticleTag::Muon);
  assert(dst_ev.particles[2].pid == dst::ParticleTag::Photon);
  assert(std::fabs(dst_ev.visible_energy_gev - 91.7f) < 0.1f);
  assert(dst_ev.to_json().find("\"visible_energy_gev\":") != std::string::npos);

  std::cout << "test_rich_vertex_dst: PASS\n";
  return 0;
}
