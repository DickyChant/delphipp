// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#include "delphipp/dst/dst_fix.hpp"
#include "delphipp/dst/energy_flow.hpp"
#include "delphipp/dst/event_assoc.hpp"
#include "delphipp/dst/short_dst.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/reco/stw.hpp"
#include "delphipp/reco/v0_finder.hpp"

int main() {
  using namespace delphipp;

  // 1. Standalone ID Straw Tube (STW) detector raw unpacker -> element reconstruction -> TER
  {
    const std::vector<std::int32_t> raw_words = {
        7, 0, 0, 0,
        static_cast<std::int32_t>((0u << 16) | (20u << 8) | 12u),
        static_cast<std::int32_t>((2u << 16) | (20u << 8) | 18u),
        static_cast<std::int32_t>((4u << 16) | (20u << 8) | 24u),
    };
    const auto hits = reco::stw::unpack_raw_blocklet(raw_words);
    assert(hits.size() == 3);
    assert(hits[0].sector == 2);
    const auto elems = reco::stw::reconstruct_elements(hits);
    assert(elems.size() == 1);
    assert(elems[0].num_layers == 3);
    TanagraStore store;
    const auto ids = reco::stw::write_ter_banks(store, elems);
    assert(ids.size() == 1);
    assert(store.ter_count() == 1);
  }

  // 2. V0 / secondary vertex finder on displaced K_S^0 -> pi+ pi- at r = 10 cm
  {
    MagneticField bf;
    bf.load_delphi_map();
    const float bz = bf.evaluate({0.0f, 0.0f, 0.0f}).bz;
    const float pt = 0.84f;
    const float inv_r_mag = bz / pt;
    const float s_approx = 9.69f;
    const float alpha = 0.25f;

    const PerigeeTrack t_pos{
        .eps = -10.0f * std::sin(alpha),
        .z0 = -s_approx / std::tan(1.50f),
        .theta = 1.50f,
        .phi0 = alpha - 0.5f * inv_r_mag * s_approx,
        .inv_r = +inv_r_mag,
    };
    const PerigeeTrack t_neg{
        .eps = +10.0f * std::sin(alpha),
        .z0 = -s_approx / std::tan(1.64f),
        .theta = 1.64f,
        .phi0 = -alpha + 0.5f * inv_r_mag * s_approx,
        .inv_r = -inv_r_mag,
    };

    const std::vector<PerigeeTrack> tracks = {t_pos, t_neg};
    const std::vector<std::int32_t> charges = {+1, -1};
    const auto v0s = reco::tracking::find_v0_candidates(tracks, charges, bz);
    assert(v0s.size() == 1);
    assert(v0s[0].dca_cm < 0.25f);
    assert(v0s[0].decay_length_cm > 8.0f && v0s[0].decay_length_cm < 12.0f);
    assert(std::fabs(v0s[0].mass_kshort_gev - 0.4976f) < 0.03f);
    assert(v0s[0].best_hypothesis == reco::tracking::V0Hypothesis::KShort);
  }

  // 3. Durham jet clustering & event shape (thrust, visible mass, y_32) on Z -> q qbar
  {
    std::vector<dst::DstParticle> ev_particles = {
        {.id = 1, .charge = +1, .px = 20.0f, .py = 1.0f, .pz = 0.5f, .energy = 20.03f},
        {.id = 2, .charge = 0,  .px = 25.0f, .py = -0.8f, .pz = -0.4f, .energy = 25.02f},
        {.id = 3, .charge = -1, .px = -22.0f, .py = -0.6f, .pz = 0.2f, .energy = 22.01f},
        {.id = 4, .charge = 0,  .px = -23.0f, .py = 0.4f, .pz = -0.3f, .energy = 23.01f},
    };
    const auto jets = dst::cluster_durham_jets(ev_particles, 0.01f, 2);
    assert(jets.size() == 2);
    const auto shape = dst::compute_event_shape(ev_particles);
    assert(shape.thrust > 0.95f);
    assert(shape.visible_mass_gev > 89.0f && shape.visible_mass_gev < 92.0f);
    assert(shape.missing_pt_gev < 0.5f);

    // 4. MC truth-to-reco association
    const auto matches = dst::associate_particles(ev_particles, ev_particles);
    assert(matches.size() == 4);
    for (const auto& m : matches) {
      assert(m.reco_idx == m.truth_idx);
      assert(m.delta_r < 1e-4f);
    }

    // 5. Post-DELANA DSTFIX / ECORR scaling & noise pruning
    dst::DstEvent dst_ev;
    dst_ev.particles = ev_particles;
    dst_ev.particles.push_back(
        dst::DstParticle{.id = 5, .charge = 0, .px = 0.02f, .py = 0.0f, .pz = 0.0f,
                         .energy = 0.02f, .pid = dst::ParticleTag::Photon});
    dst::apply_dst_fixes(dst_ev, dst::DstFixConfig{.ecal_scale = 1.05f, .tpc_z_shift_cm = 0.0232f,
                                                   .min_neutral_energy_gev = 0.10f});
    assert(dst_ev.particles.size() == 4);  // noise cluster pruned
    assert(std::fabs(dst_ev.particles[0].impact_z_cm - 0.0232f) < 1e-5f);
    assert(std::fabs(dst_ev.particles[1].energy - 25.02f * 1.05f) < 1e-2f);
    assert(dst_ev.visible_energy_gev > 90.0f);
  }

  std::cout << "test_physics_dst: PASS\n";
  return 0;
}
