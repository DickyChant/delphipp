// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// DELPHI++ Vertical Slice v0.2 end-to-end pure C++20 event pipeline driver.
// Exercises db/, core/ (MagneticField, PerigeeTrack, TanagraStore), io/
// (FzExchangeWriter/Reader), and reco/ (TOF, MUS, MUF, MUB), emitting
// bank_diff.py-compatible JSON output.

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "delphipp/db/condition_db.hpp"
#include "delphipp/helix.hpp"
#include "delphipp/io/fz_exchange.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/reco/mub.hpp"
#include "delphipp/reco/muf.hpp"
#include "delphipp/reco/mus.hpp"
#include "delphipp/reco/tof.hpp"
#include "delphipp/tanagra_store.hpp"

int main(int argc, char** argv) {
  std::string out_json_path;
  std::int32_t run_number = 47014;
  std::int32_t event_number = 1;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--out" && i + 1 < argc) {
      out_json_path = argv[++i];
    } else if (arg == "--run" && i + 1 < argc) {
      run_number = std::stoi(argv[++i]);
    } else if (arg == "--event" && i + 1 < argc) {
      event_number = std::stoi(argv[++i]);
    }
  }

  // 1. Condition DB (94D TPC ABST0 correction)
  delphipp::db::ConditionDatabase cond_db;
  cond_db.add_record(delphipp::db::ConditionRecord{
      .table = "CALB",
      .path = "GLOB$TPC*.GLOB$ENP0",
      .is_word_update = false,
      .valid_from = {931101, 0},
      .entered_on = {940601, 120000},
      .field = "SLOW",
      .words = {{1, 6.72f}, {14, 1.147f}},
  });
  std::istringstream asc(
      "*DBRECORD\nCALB GLOB$TPC*.GLOB$ENP0.W\n 940101,000001,990101,120000\n"
      "*DBFIELD SLOW\n 3,SLOW,14., 1.10208\n");
  cond_db.load_ascii(asc);
  const auto abst0 = cond_db.resolve_word("CALB", "GLOB$TPC*.GLOB$ENP0", "SLOW",
                                          14, {940526, 120000});

  // 2. Magnetic field & Z -> mu+ mu- (45.6 GeV barrel muon) helix propagation
  delphipp::MagneticField bfield;
  const auto b0 = bfield.evaluate({0.0f, 0.0f, 0.0f});
  const float p_mu = 45.6f;
  const float inv_r = b0.bz / p_mu;
  delphipp::PerigeeTrack mu_track{
      .eps = 0.0f,
      .z0 = 0.0f,
      .theta = 1.45f,
      .phi0 = 0.35f,
      .inv_r = inv_r,
  };

  const auto at_tof = delphipp::propagate_to_cylinder(mu_track, delphipp::reco::tof::kTofRadiusCm);
  const auto at_mub = delphipp::propagate_to_cylinder(mu_track, 445.0f);

  // 3. Reconstruct into unified TanagraStore
  delphipp::TanagraStore store;

  // TOF hit at extrapolated barrel crossing
  if (at_tof.has_value()) {
    const float r = std::hypot(at_tof->xyz[0], at_tof->xyz[1]);
    const float phi = std::atan2(at_tof->xyz[1], at_tof->xyz[0]);
    std::vector<delphipp::reco::tof::TofCounterHit> tof_hits = {
        {.counter_num = 12,
         .r_cm = r,
         .r_phi_cm = r * phi,
         .z_cm = at_tof->xyz[2],
         .pm_ext = 15.2f,
         .pm_int = 14.8f},
    };
    delphipp::reco::tof::create_td_banks(store, tof_hits);
  }

  // MUB triplet hit
  delphipp::reco::mub::MubChamberCalib mub_ch{};
  mub_ch.tzero_ns = {100.0f, 100.0f, 100.0f};
  mub_ch.delay_poly = {0.0f, 0.5f, 0.0f, 0.0f, 0.0f};
  const auto mub_hit =
      delphipp::reco::mub::reconstruct_triplet(1, {207.3f, 500.0f, 500.0f}, mub_ch);
  if (mub_hit.has_value()) {
    const auto td_mub = store.add_tdr(delphipp::TdrBank{
        .detector_id = delphipp::reco::mub::kDetectorId,
        .submod_id = 1,
        .descriptor = 121,
        .num_points = 1,
        .words_per_point = 4,
        .data = {mub_hit->z_local_cm, mub_hit->drift_cm,
                 mub_hit->anode_time_ns, mub_hit->time_sum_ns},
    });
    if (at_mub.has_value()) {
      store.add_ter(delphipp::TerBank{
          .detector_id = delphipp::reco::mub::kDetectorId,
          .submod_id = 1,
          .descriptor = 125,
          .fwords = {at_mub->xyz[0], at_mub->xyz[1], at_mub->xyz[2],
                     abst0.value_or(0.0f)},
          .td_refs = {td_mub},
      });
    }
  }

  // 4. Round-trip event summary through FZ Exchange binary stream
  std::stringstream fz_stream(std::ios::in | std::ios::out | std::ios::binary);
  {
    delphipp::io::FzExchangeWriter writer(fz_stream, 64);
    writer.write_record(delphipp::io::FzLogicalRecord{
        .type = delphipp::io::FzRecordType::EventNormal,
        .user_header = {run_number, event_number},
        .payload = {static_cast<std::int32_t>(store.tdr_count()),
                    static_cast<std::int32_t>(store.ter_count())},
    });
    writer.finish();
  }
  fz_stream.seekg(0);
  delphipp::io::FzExchangeReader reader(fz_stream);
  const auto fz_rec = reader.next_record();
  if (!fz_rec.has_value() || fz_rec->user_header[0] != run_number) {
    std::cerr << "delphipp_reco: FZ round-trip verification failed\n";
    return 1;
  }

  // 5. Emit bank_diff.py-compatible JSON
  const std::string json = store.to_bank_diff_json(run_number, event_number);
  if (!out_json_path.empty()) {
    std::ofstream ofs(out_json_path);
    ofs << json << '\n';
  } else {
    std::cout << json << '\n';
  }
  return 0;
}
