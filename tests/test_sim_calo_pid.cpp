// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "delphipp/io/raw_blocklet.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/pipeline/event_pipeline.hpp"
#include "delphipp/reco/emf.hpp"
#include "delphipp/reco/had.hpp"
#include "delphipp/reco/hof.hpp"
#include "delphipp/reco/hpc.hpp"
#include "delphipp/reco/mub.hpp"
#include "delphipp/reco/muf.hpp"
#include "delphipp/reco/mus.hpp"
#include "delphipp/reco/rich.hpp"
#include "delphipp/reco/sat.hpp"
#include "delphipp/reco/stc.hpp"
#include "delphipp/reco/tof.hpp"
#include "delphipp/sim/calo_pid_sim.hpp"
#include "delphipp/zebra_store.hpp"

namespace {

void check_cond(bool cond, const char* expr, int line) {
  if (!cond) {
    std::cerr << "test_sim_calo_pid FAILED at line " << line << ": " << expr << '\n';
    std::abort();
  }
}

#define CHECK(expr) check_cond((expr), #expr, __LINE__)

}  // namespace

int main() {
  using namespace delphipp;

  MagneticField bfield;
  bfield.load_delphi_map();
  pipeline::EventReconstructionPipeline pipe(bfield, 1.10208f);

  // 1. Multi-particle event covering barrel electron, endcap electron, muon,
  //    neutral hadron, and forward Bhabha e+e- pair
  const std::vector<sim::SimShowerParticle> particles = {
      // Barrel electron (|pdg|=11, theta ~ 1.46 rad -> HPC, TOF, BRICH)
      {.pdg = 11, .charge = -1, .px = 15.0f, .py = 10.0f, .pz = 2.0f, .energy = 18.14f},
      // Endcap positron (|pdg|=11, theta ~ 0.28 rad -> EMF, FRICH, HOF)
      {.pdg = -11, .charge = 1, .px = 6.0f, .py = -2.0f, .pz = 22.0f, .energy = 22.89f},
      // Muon in barrel-forward surround overlap (|pdg|=13, theta ~ 0.785 rad -> MUB, MUF, MUS, TOF, HOF, BRICH, HAC MIP)
      {.pdg = 13, .charge = -1, .px = 12.0f, .py = 12.0f, .pz = 17.0f, .energy = 23.98f},
      // Neutral hadron (pdg=2112, theta ~ 1.38 rad, phi ~ 2.88 rad -> HAC)
      {.pdg = 2112, .charge = 0, .px = -15.0f, .py = 4.0f, .pz = 3.0f, .energy = 15.84f},
      // Forward Bhabha pair (|pdg|=11, back-to-back at theta ~ 0.09 rad & pi - 0.09 rad -> STC, SAT, LUM)
      {.pdg = 11, .charge = -1, .px = 3.6f, .py = 1.8f, .pz = 45.4f, .energy = 45.6f},
      {.pdg = -11, .charge = 1, .px = -3.6f, .py = -1.8f, .pz = -45.4f, .energy = 45.6f},
  };

  const auto fz_ev = sim::digitize_calo_muon_rich(-47014, 17, particles, 42ULL);
  CHECK(fz_ev.run_number == -47014);
  CHECK(fz_ev.event_number == 17);
  CHECK(fz_ev.has_header_event_number);
  CHECK(fz_ev.top_bank == "TOP ");

  // Determinism check for identical seed
  const auto fz_ev_repeat = sim::digitize_calo_muon_rich(-47014, 17, particles, 42ULL);
  CHECK(fz_ev.banks.size() == fz_ev_repeat.banks.size());
  for (std::size_t i = 0; i < fz_ev.banks.size(); ++i) {
    CHECK(fz_ev.banks[i].idh == fz_ev_repeat.banks[i].idh);
    CHECK(fz_ev.banks[i].iwords == fz_ev_repeat.banks[i].iwords);
  }

  // Pass through unified EventReconstructionPipeline
  const auto out = pipe.reconstruct_event(fz_ev);
  CHECK(out.run == -47014);
  CHECK(out.event == 17);
  CHECK(out.top_bank == "TOP ");
  CHECK(out.hof_hits >= 1);
  CHECK(out.dst_event.is_bhabha);
  CHECK(std::isfinite(out.dst_event.visible_energy_gev));
  CHECK(out.dst_event.visible_energy_gev > 100.0f);
  CHECK(!out.dst_event.particles.empty());

  for (const char* det : {"HPC", "EMF", "HAC", "STC", "SAT", "LUM", "TOF",
                          "HOF", "MUB", "MUF", "MUS", "RIB", "RIF"}) {
    CHECK(out.detector_blocklet_counts.contains(det));
    CHECK(out.detector_blocklet_counts.at(det) >= 1);
  }

  // Verify TANAGRA reconstructed TER/TDR banks exist for each subdetector
  auto has_ter_with_desc = [&](std::int32_t det_id, std::int32_t desc) {
    for (const auto id : out.tanagra.list_ter(det_id)) {
      const auto* tb = out.tanagra.get_ter(id);
      if (tb != nullptr && tb->descriptor == desc) {
        return true;
      }
    }
    return false;
  };
  CHECK(has_ter_with_desc(reco::hpc::kDetectorId, 125));
  CHECK(has_ter_with_desc(reco::emf::kDetectorId, 124));
  CHECK(has_ter_with_desc(reco::had::kDetectorId, 126));
  CHECK(has_ter_with_desc(reco::stc::kDetectorId, 125));
  CHECK(has_ter_with_desc(reco::sat::kDetectorIdSat, 125));
  CHECK(has_ter_with_desc(reco::sat::kDetectorIdVsat, 125));
  CHECK(has_ter_with_desc(reco::mub::kDetectorId, 125));
  CHECK(has_ter_with_desc(reco::muf::kDetectorId, 125));
  CHECK(has_ter_with_desc(reco::mus::kDetectorId, 125));

  bool has_tof_tdr = false;
  for (const auto id : out.tanagra.list_tdr(reco::tof::kDetectorId)) {
    const auto* td = out.tanagra.get_tdr(id);
    if (td != nullptr && td->words_per_point == 9) {
      has_tof_tdr = true;
      break;
    }
  }
  CHECK(has_tof_tdr);

  bool has_hof_tdr = false;
  for (const auto id : out.tanagra.list_tdr(reco::hof::kDetectorId)) {
    const auto* td = out.tanagra.get_tdr(id);
    if (td != nullptr && td->descriptor == 290 && td->words_per_point == 8) {
      has_hof_tdr = true;
      break;
    }
  }
  CHECK(has_hof_tdr);

  // Direct ZEBRA traversal & unpacker/clusterer verification
  OwnedZebraStore store(65536);
  auto view = store.view();
  const auto heads = fz_ev.populate_zebra_store(view);
  CHECK(heads.contains("TOP "));
  const std::int32_t lraw = view.bank_link(heads.at("TOP "), 1);
  CHECK(lraw > 0);

  auto read_part_blocklets = [&](std::int32_t raw_link_idx) {
    const std::int32_t lpmb = view.bank_link(lraw, raw_link_idx);
    CHECK(lpmb > 0);
    const std::int32_t lrawd = view.bank_link(lpmb, 5);
    CHECK(lrawd > 0);
    return io::read_raw_blocklet_chain(view, lrawd);
  };

  // HPC: 9 longitudinal layers -> 1 shower
  {
    const auto blks = read_part_blocklets(11);
    CHECK(blks.size() == 1);
    const auto hits = reco::hpc::unpack_raw_blocklet(blks[0].data);
    CHECK(hits.size() >= 9);
    const auto showers = reco::hpc::reconstruct_showers(hits);
    CHECK(showers.size() == 1);
    CHECK(showers[0].energy_gev > 15.0f && showers[0].energy_gev < 21.0f);
  }

  // EMF: lead-glass 3x3 cluster -> 1 cluster
  {
    const auto blks = read_part_blocklets(19);
    CHECK(blks.size() == 1);
    const auto blocks = reco::emf::unpack_raw_blocklet(blks[0].data);
    CHECK(!blocks.empty());
    const auto clusters = reco::emf::find_clusters(blocks);
    CHECK(clusters.size() == 1);
    CHECK(clusters[0].energy_gev > 18.0f && clusters[0].energy_gev < 26.0f);
  }

  // HAC: 0x5A5A marker + 4 layers -> hadron + muon MIP clusters
  {
    const auto blks = read_part_blocklets(13);
    CHECK(blks.size() == 1);
    CHECK(blks[0].data.size() >= 6);
    CHECK((static_cast<std::uint32_t>(blks[0].data[2]) & 0xFFFFu) == 0x5A5Au);
    const auto towers = reco::had::unpack_raw_blocklet(blks[0].data);
    CHECK(towers.size() >= 8);
    const auto clusters = reco::had::find_clusters(towers, 0.25f);
    CHECK(clusters.size() >= 2);
    CHECK(clusters[0].energy_gev > 12.0f);
  }

  // STC (STRA) + SAT + LUM: forward Bhabha arms
  {
    CHECK(heads.contains("STRA"));
    const auto stc_blks = io::read_raw_blocklet_chain(view, heads.at("STRA"));
    CHECK(stc_blks.size() == 1);
    const auto stc_towers = reco::stc::unpack_raw_blocklet(stc_blks[0].data);
    CHECK(stc_towers.size() == 2);
    const auto stc_clusters = reco::stc::find_clusters(stc_towers);
    CHECK(stc_clusters.size() == 2);
    CHECK(reco::stc::is_bhabha_candidate(stc_clusters));

    const auto sat_blks = read_part_blocklets(17);
    CHECK(sat_blks.size() == 1);
    const auto sat_hits = reco::sat::unpack_raw_blocklet(sat_blks[0].data);
    CHECK(sat_hits.size() == 2);
    const auto sat_clusters = reco::sat::reconstruct_clusters(sat_hits);
    CHECK(sat_clusters.size() == 2);

    const auto lum_blks = read_part_blocklets(22);
    CHECK(lum_blks.size() == 1);
    const auto lum_hits = reco::sat::unpack_raw_blocklet(lum_blks[0].data);
    CHECK(lum_hits.size() == 2);
  }

  // TOF & HOF
  {
    const auto tof_blks = read_part_blocklets(12);
    CHECK(tof_blks.size() == 1);
    const auto tof_hits = reco::tof::unpack_raw_blocklet(tof_blks[0].data);
    CHECK(tof_hits.size() == 2);

    const std::int32_t lpmb_hof = view.bank_link(lraw, 20);
    CHECK(lpmb_hof > 0);
    const std::int32_t lrawd_hof = view.bank_link(lpmb_hof, 5);
    CHECK(lrawd_hof > 0);
    reco::hof::HofCalibration hcal{};
    reco::hof::HofStats hstats{};
    hstats.begin_event();
    const auto hof_td = reco::hof::unpack_td(view, lrawd_hof, 0, hcal, hstats);
    CHECK(hof_td.counter_hits >= 1);
  }

  // MUB, MUF, MUS (EDRA strip+tube pairs -> MusTeBank)
  {
    const auto mub_blks = read_part_blocklets(14);
    CHECK(mub_blks.size() == 1);
    const auto mub_hits = reco::mub::unpack_raw_blocklet(mub_blks[0].data);
    CHECK(mub_hits.size() == 4);

    const auto muf_blks = read_part_blocklets(21);
    CHECK(muf_blks.size() == 1);
    const auto muf_trips = reco::muf::unpack_raw_blocklet(muf_blks[0].data);
    CHECK(muf_trips.size() == 2);

    CHECK(heads.contains("EDRA"));
    const auto edra_blks = io::read_raw_blocklet_chain(view, heads.at("EDRA"));
    CHECK(edra_blks.size() == 1);
    const auto mus_raw = reco::mus::unpack_raw_blocklet(edra_blks[0].data);
    CHECK(mus_raw.size() == 4);
    const auto mus_clus = reco::mus::find_clusters(mus_raw);
    CHECK(mus_clus.size() == 4);
    reco::mus::MusGeometry mgeom{};
    for (std::size_t d = 0; d < reco::mus::kNumDetectors; ++d) {
      for (std::size_t l = 0; l < reco::mus::kNumLayers; ++l) {
        mgeom.pitch_ts1[d][l] = {5.0f, 5.0f};
      }
      mgeom.local_to_delphi[d][0] = {470.0f, 0.0f, 180.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f};
      mgeom.local_to_delphi[d][1] = {470.0f, 0.0f, 180.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f};
    }
    const auto mus_sp = reco::mus::combine_clusters(mus_clus, mgeom);
    CHECK(mus_sp.size() == 2);
    const auto mus_te = reco::mus::build_te(mus_sp, mgeom);
    CHECK(mus_te.size() == 2);
  }

  // BRICH (RIB) & FRICH (RIF) Cherenkov photons
  {
    const auto rib_blks = read_part_blocklets(9);
    CHECK(rib_blks.size() == 1);
    const auto rib_ph = reco::rich::unpack_raw_blocklet(rib_blks[0].data, reco::rich::kDetectorIdBrich);
    CHECK(!rib_ph.empty());
    for (const auto& ph : rib_ph) {
      CHECK(ph.theta_c_rad > 0.60f && ph.theta_c_rad <= 0.66f);
    }

    const auto rif_blks = read_part_blocklets(16);
    CHECK(rif_blks.size() == 1);
    const auto rif_ph = reco::rich::unpack_raw_blocklet(rif_blks[0].data, reco::rich::kDetectorIdFrich);
    CHECK(!rif_ph.empty());
    for (const auto& ph : rif_ph) {
      CHECK(ph.theta_c_rad > 0.60f && ph.theta_c_rad <= 0.66f);
    }
  }

  // 2. Empty particle event edge case
  {
    const auto empty_ev = sim::digitize_calo_muon_rich(47014, 1, {}, 1ULL);
    const auto empty_out = pipe.reconstruct_event(empty_ev);
    CHECK(empty_out.raw_blocklets == 0);
    CHECK(empty_out.hof_hits == 0);
    CHECK(std::isfinite(empty_out.dst_event.visible_energy_gev));
    CHECK(empty_out.dst_event.visible_energy_gev == 0.0f);
  }

  std::cout << "test_sim_calo_pid: PASS\n";
  return 0;
}
