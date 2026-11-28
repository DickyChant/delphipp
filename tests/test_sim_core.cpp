// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "delphipp/io/raw_blocklet.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/pipeline/event_pipeline.hpp"
#include "delphipp/reco/fwd_chambers.hpp"
#include "delphipp/reco/id.hpp"
#include "delphipp/reco/od.hpp"
#include "delphipp/reco/stw.hpp"
#include "delphipp/reco/tpc.hpp"
#include "delphipp/reco/vd.hpp"
#include "delphipp/sim/sim_core.hpp"
#include "delphipp/sim/tracker_sim.hpp"
#include "delphipp/zebra_store.hpp"

int main() {
  using namespace delphipp;

  // 1. Deterministic LEP1 beamspot smearing
  {
    sim::PrimaryVertexSmearer smear_a(123456ULL);
    sim::PrimaryVertexSmearer smear_b(123456ULL);
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double sq_x = 0.0, sq_y = 0.0, sq_z = 0.0;
    constexpr int kSamples = 2000;
    for (int i = 0; i < kSamples; ++i) {
      const Vec3 va = smear_a.sample();
      const Vec3 vb = smear_b.sample();
      assert(va == vb);
      for (float c : va) {
        assert(std::isfinite(c));
      }
      sum_x += va[0];
      sum_y += va[1];
      sum_z += va[2];
      sq_x += static_cast<double>(va[0]) * va[0];
      sq_y += static_cast<double>(va[1]) * va[1];
      sq_z += static_cast<double>(va[2]) * va[2];
    }
    const double mean_x = sum_x / kSamples;
    const double mean_y = sum_y / kSamples;
    const double mean_z = sum_z / kSamples;
    const double rms_x = std::sqrt(sq_x / kSamples - mean_x * mean_x);
    const double rms_y = std::sqrt(sq_y / kSamples - mean_y * mean_y);
    const double rms_z = std::sqrt(sq_z / kSamples - mean_z * mean_z);
    assert(std::fabs(mean_x) < 0.003);
    assert(std::fabs(mean_y) < 0.0003);
    assert(std::fabs(mean_z) < 0.2);
    assert(std::fabs(rms_x - 0.015) < 0.003);
    assert(std::fabs(rms_y - 0.001) < 0.0003);
    assert(std::fabs(rms_z - 1.0) < 0.2);
  }

  MagneticField bfield;
  bfield.load_delphi_map();

  // 2. B-field stepping across DELPHI tracker layers
  {
    const sim::SimParticle neutral{
        .pdg = 22,
        .charge = 0,
        .px = 5.0f,
        .py = 3.0f,
        .pz = 1.0f,
        .energy = 5.916f,
        .mass = 0.0f,
        .vertex_cm = {0.0f, 0.0f, 0.0f},
    };
    assert(sim::step_particle_in_bfield(neutral, bfield, 1ULL).empty());

    const sim::SimParticle barrel_pi{
        .pdg = 211,
        .charge = +1,
        .px = 5.0f * std::cos(0.52f),
        .py = 5.0f * std::sin(0.52f),
        .pz = 1.1f,
        .energy = 5.121f,
        .mass = 0.13957f,
        .vertex_cm = {0.01f, -0.001f, 0.5f},
    };
    const auto hits_a = sim::step_particle_in_bfield(barrel_pi, bfield, 77ULL);
    const auto hits_b = sim::step_particle_in_bfield(barrel_pi, bfield, 77ULL);
    assert(hits_a.size() == hits_b.size());
    std::set<std::string> seen_barrel;
    for (std::size_t i = 0; i < hits_a.size(); ++i) {
      assert(hits_a[i].det_name == hits_b[i].det_name);
      assert(hits_a[i].r_cm == hits_b[i].r_cm);
      assert(hits_a[i].phi_rad == hits_b[i].phi_rad);
      assert(hits_a[i].z_cm == hits_b[i].z_cm);
      assert(std::isfinite(hits_a[i].path_x0) && hits_a[i].path_x0 > 0.0f);
      assert(std::isfinite(hits_a[i].dedx_mev) && hits_a[i].dedx_mev > 0.0f);
      seen_barrel.insert(hits_a[i].det_name);
    }
    for (const char* det : {"VD", "ID", "STW", "TPC", "OD"}) {
      assert(seen_barrel.contains(det));
    }

    const float pt_fwd = 3.0f;
    const float pz_fwd = pt_fwd / std::tan(0.35f);
    const sim::SimParticle fwd_mu{
        .pdg = -13,
        .charge = +1,
        .px = pt_fwd * std::cos(2.10f),
        .py = pt_fwd * std::sin(2.10f),
        .pz = pz_fwd,
        .energy = std::hypot(std::hypot(pt_fwd, pz_fwd), 0.10566f),
        .mass = 0.10566f,
        .vertex_cm = {0.0f, 0.0f, 0.0f},
    };
    const auto fwd_hits = sim::step_particle_in_bfield(fwd_mu, bfield, 88ULL);
    std::set<std::string> seen_fwd;
    for (const auto& h : fwd_hits) {
      seen_fwd.insert(h.det_name);
    }
    assert(seen_fwd.contains("FCA"));
    assert(seen_fwd.contains("FCB"));
  }

  // 3. Full SimEvent -> digitize_trackers -> 6 tracker unpackers + Pipeline -> ShortDST
  sim::PrimaryVertexSmearer vtx_smear(2026ULL);
  const Vec3 pv = vtx_smear.sample();
  sim::SimEvent sim_ev{
      .run = -47014,
      .event = 42,
      .primary_vertex_cm = pv,
      .particles = {
          sim::SimParticle{
              .pdg = 211,
              .charge = +1,
              .px = 5.0f * std::cos(0.52f),
              .py = 5.0f * std::sin(0.52f),
              .pz = 1.10f,
              .energy = 5.121f,
              .mass = 0.13957f,
              .vertex_cm = pv,
          },
          sim::SimParticle{
              .pdg = -211,
              .charge = -1,
              .px = 6.0f * std::cos(3.66f),
              .py = 6.0f * std::sin(3.66f),
              .pz = -1.20f,
              .energy = 6.120f,
              .mass = 0.13957f,
              .vertex_cm = pv,
          },
          sim::SimParticle{
              .pdg = -13,
              .charge = +1,
              .px = 3.0f * std::cos(2.10f),
              .py = 3.0f * std::sin(2.10f),
              .pz = 3.0f / std::tan(0.35f),
              .energy = 8.75f,
              .mass = 0.10566f,
              .vertex_cm = pv,
          },
      },
  };

  const io::ParsedFzEvent fz_ev = sim::digitize_trackers(sim_ev, bfield, 9001ULL);
  assert(fz_ev.top_bank == "TOP ");

  // Verify 1-to-1 decoding through all 6 tracker unpackers via relocated ZEBRA store
  {
    OwnedZebraStore zstore(65536);
    auto zview = zstore.view();
    const auto heads = fz_ev.populate_zebra_store(zview);
    assert(heads.contains("TOP "));
    const std::int32_t lraw = zview.bank_link(heads.at("TOP "), 1);
    assert(lraw > 0);

    auto read_pmb_data = [&](std::int32_t link_idx) {
      const std::int32_t lpmb = zview.bank_link(lraw, link_idx);
      assert(lpmb > 0);
      const std::int32_t lrawd = zview.bank_link(lpmb, 5);
      assert(lrawd > 0);
      const auto chain = io::read_raw_blocklet_chain(zview, lrawd);
      assert(chain.size() == 1);
      return chain.front().data;
    };

    const auto vd_hits = reco::vd::unpack_raw_blocklet(read_pmb_data(6));
    assert(!reco::vd::find_clusters(vd_hits).empty());

    const auto id_hits = reco::id::unpack_raw_blocklet(read_pmb_data(7));
    assert(!reco::id::reconstruct_segments(id_hits).empty());

    const auto tpc_hits = reco::tpc::unpack_raw_blocklet(read_pmb_data(8));
    assert(!reco::tpc::reconstruct_space_points(tpc_hits, {}).empty());

    const auto od_hits = reco::od::unpack_raw_blocklet(read_pmb_data(10));
    assert(!reco::od::reconstruct_planks(od_hits).empty());

    const auto fca_hits = reco::fwd::unpack_raw_blocklet(
        read_pmb_data(15), reco::fwd::kDetectorIdFca);
    assert(reco::fwd::fit_segment(fca_hits).has_value());

    const auto fcb_hits = reco::fwd::unpack_raw_blocklet(
        read_pmb_data(18), reco::fwd::kDetectorIdFcb);
    assert(reco::fwd::fit_segment(fcb_hits).has_value());

    const auto* stwr = fz_ev.find_first("STWR");
    assert(stwr != nullptr);
    const auto stw_hits = reco::stw::unpack_raw_blocklet(stwr->iwords);
    assert(!reco::stw::reconstruct_elements(stw_hits).empty());
  }

  // Closed-loop reconstruction pipeline -> ShortDST
  pipeline::EventReconstructionPipeline pipe(bfield, 1.10208f);
  const auto reco_out = pipe.reconstruct_event(fz_ev);
  assert(reco_out.run == -47014);
  assert(reco_out.event == 42);
  for (const char* det : {"VD", "ID", "TPC", "OD", "FCA", "FCB"}) {
    assert(reco_out.detector_blocklet_counts.contains(det));
    assert(reco_out.detector_blocklet_counts.at(det) >= 1);
  }
  assert(reco_out.kalman_refits >= 1);
  assert(reco_out.tracks_after_ambiguity >= 1);
  assert(reco_out.tanagra.tkr_count() >= 1);
  assert(!reco_out.dst_event.particles.empty());
  assert(std::isfinite(reco_out.dst_event.visible_energy_gev));
  assert(reco_out.dst_event.visible_energy_gev > 1.0f);
  for (const auto& p : reco_out.dst_event.particles) {
    assert(std::isfinite(p.px) && std::isfinite(p.py) && std::isfinite(p.pz));
    assert(std::isfinite(p.energy));
    assert(std::isfinite(p.impact_rphi_cm) && std::isfinite(p.impact_z_cm));
  }

  // 4. SIM_ truth record builder
  {
    const auto truth_ev = sim::build_sim_truth_record(sim_ev);
    assert(truth_ev.top_bank == "SIM_");
    assert(truth_ev.find_first("SIM_") != nullptr);
    assert(truth_ev.find_first("SP  ") != nullptr);
    assert(truth_ev.find_all("ST  ").size() == 3);
    const auto truth_out = pipe.reconstruct_event(truth_ev);
    assert(truth_out.dst_event.particles.size() == 3);
    assert(std::isfinite(truth_out.dst_event.visible_energy_gev));
    assert(truth_out.dst_event.visible_energy_gev > 15.0f);
  }

  std::cout << "test_sim_core: PASS\n";
  return 0;
}
