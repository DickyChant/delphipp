// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// DELPHI++ pure-C++20 reconstruction driver.
// Supports both:
//   1. Self-contained synthetic Z -> mu+ mu- event mode (when --input is omitted).
//   2. Real CERNLIB ZEBRA FZ Exchange file mode (--input <file.fadsim|.sdst|.fz>):
//      runs EventReconstructionPipeline across parsed events and exports
//      bank_diff JSON, structured summary JSON (--summary-json), CSV, and HepMC3.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "delphipp/db/condition_db.hpp"
#include "delphipp/dst/dst_writer.hpp"
#include "delphipp/dst/short_dst.hpp"
#include "delphipp/helix.hpp"
#include "delphipp/io/fz_exchange.hpp"
#include "delphipp/io/fz_zebra_parser.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/pipeline/event_pipeline.hpp"
#include "delphipp/reco/mub.hpp"
#include "delphipp/reco/tof.hpp"
#include "delphipp/tanagra_store.hpp"

namespace {

struct CliOptions {
  std::string input_fz_path;
  std::string out_json_path;
  std::string summary_json_path;
  std::string csv_out_path;
  std::string hepmc_out_path;
  std::int32_t run_number = 47014;
  std::int32_t event_number = 1;
  std::size_t max_events = 0;
  bool verbose = false;
};

CliOptions parse_cli(int argc, char** argv) {
  CliOptions opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--input" && i + 1 < argc) {
      opt.input_fz_path = argv[++i];
    } else if (arg == "--out" && i + 1 < argc) {
      opt.out_json_path = argv[++i];
    } else if (arg == "--summary-json" && i + 1 < argc) {
      opt.summary_json_path = argv[++i];
    } else if (arg == "--csv-out" && i + 1 < argc) {
      opt.csv_out_path = argv[++i];
    } else if (arg == "--hepmc-out" && i + 1 < argc) {
      opt.hepmc_out_path = argv[++i];
    } else if (arg == "--run" && i + 1 < argc) {
      opt.run_number = std::stoi(argv[++i]);
    } else if (arg == "--event" && i + 1 < argc) {
      opt.event_number = std::stoi(argv[++i]);
    } else if (arg == "--max-events" && i + 1 < argc) {
      opt.max_events = static_cast<std::size_t>(std::stoul(argv[++i]));
    } else if (arg == "--verbose" || arg == "-v") {
      opt.verbose = true;
    }
  }
  return opt;
}

int run_synthetic_event(const CliOptions& opt,
                        float abst0_val,
                        const delphipp::MagneticField& bfield) {
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

  const auto at_tof =
      delphipp::propagate_to_cylinder(mu_track, delphipp::reco::tof::kTofRadiusCm);
  const auto at_mub = delphipp::propagate_to_cylinder(mu_track, 445.0f);

  delphipp::TanagraStore store;
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
          .fwords = {at_mub->xyz[0], at_mub->xyz[1], at_mub->xyz[2], abst0_val},
          .td_refs = {td_mub},
      });
    }
  }

  std::stringstream fz_stream(std::ios::in | std::ios::out | std::ios::binary);
  {
    delphipp::io::FzExchangeWriter writer(fz_stream, 64);
    writer.write_record(delphipp::io::FzLogicalRecord{
        .type = delphipp::io::FzRecordType::EventNormal,
        .user_header = {opt.run_number, opt.event_number},
        .payload = {static_cast<std::int32_t>(store.tdr_count()),
                    static_cast<std::int32_t>(store.ter_count())},
    });
    writer.finish();
  }
  fz_stream.seekg(0);
  delphipp::io::FzExchangeReader reader(fz_stream);
  const auto fz_rec = reader.next_record();
  if (!fz_rec.has_value() || fz_rec->user_header[0] != opt.run_number) {
    std::cerr << "delphipp_reco: FZ round-trip verification failed\n";
    return 1;
  }

  std::vector<delphipp::dst::TrackPidInput> tracks = {
      {.track = mu_track, .charge = +1, .ecal_energy_gev = 0.3f, .hcal_energy_gev = 2.1f, .has_muon_hit = true},
  };
  const auto dst_ev = delphipp::dst::build_short_dst(
      opt.run_number, opt.event_number, b0.bz, tracks, {});

  if (!opt.csv_out_path.empty()) {
    std::ofstream csv_ofs(opt.csv_out_path);
    delphipp::dst::write_csv_header(csv_ofs);
    delphipp::dst::write_csv_event(csv_ofs, dst_ev);
  }
  if (!opt.hepmc_out_path.empty()) {
    std::ofstream hm_ofs(opt.hepmc_out_path);
    delphipp::dst::write_hepmc3_header(hm_ofs);
    delphipp::dst::write_hepmc3_event(hm_ofs, dst_ev);
    delphipp::dst::write_hepmc3_footer(hm_ofs);
  }

  const std::string json = store.to_bank_diff_json(opt.run_number, opt.event_number);
  if (opt.verbose) {
    const auto b_yoke = bfield.evaluate({350.0f, 0.0f, 100.0f});
    std::cout << "=== DELPHI++ Reco Summary (run=" << opt.run_number
              << ", event=" << opt.event_number << ") ===\n"
              << "  CondDB 94D TPC ABST0   : " << abst0_val << " time bins\n"
              << "  UFIELD center Bz       : " << b0.bz_tesla() << " T\n"
              << "  UFIELD yoke(350,0,100) : Bz=" << b_yoke.bz_tesla() << " T\n"
              << "  TanagraStore banks     : " << store.tdr_count() << " TDR, "
              << store.ter_count() << " TER\n"
              << "  FZ Exchange round-trip : OK (header run=" << fz_rec->user_header[0]
              << ", evt=" << fz_rec->user_header[1] << ")\n"
              << "  Bank JSON              : " << json << '\n';
  }
  if (!opt.out_json_path.empty()) {
    std::ofstream ofs(opt.out_json_path);
    ofs << json << '\n';
  } else if (!opt.verbose) {
    std::cout << json << '\n';
  }
  return 0;
}

int run_fz_file_reco(const CliOptions& opt,
                     float abst0_val,
                     const delphipp::MagneticField& bfield) {
  std::ifstream ifs(opt.input_fz_path, std::ios::binary);
  if (!ifs) {
    std::cerr << "delphipp_reco: cannot open input file: " << opt.input_fz_path << '\n';
    return 1;
  }

  const auto parsed = delphipp::io::parse_fz_stream(ifs, opt.max_events);
  if (!parsed.detected_zebra_fz) {
    std::cerr << "delphipp_reco: file is not a valid CERNLIB ZEBRA FZ stream: "
              << opt.input_fz_path << '\n';
    return 1;
  }
  if (parsed.events.empty()) {
    std::cerr << "delphipp_reco: no events decoded from FZ stream: "
              << opt.input_fz_path << '\n';
    return 1;
  }

  const auto b0 = bfield.evaluate({0.0f, 0.0f, 0.0f});
  delphipp::pipeline::EventReconstructionPipeline pipeline(bfield, abst0_val);

  std::ofstream csv_ofs;
  if (!opt.csv_out_path.empty()) {
    csv_ofs.open(opt.csv_out_path);
    delphipp::dst::write_csv_header(csv_ofs);
  }
  std::ofstream hm_ofs;
  if (!opt.hepmc_out_path.empty()) {
    hm_ofs.open(opt.hepmc_out_path);
    delphipp::dst::write_hepmc3_header(hm_ofs);
  }

  if (opt.verbose) {
    std::cout << "=== DELPHI++ Real FZ Reconstruction ===\n"
              << "  Input file        : " << opt.input_fz_path
              << " (" << parsed.file_bytes << " bytes)\n"
              << "  Endianness        : " << (parsed.is_big_endian ? "big" : "little")
              << ", phys_words=" << parsed.physical_record_words
              << ", steering=" << parsed.steering_blocks
              << ", fast=" << parsed.fast_blocks << '\n'
              << "  Parsed events     : " << parsed.events.size()
              << ", total ZEBRA banks=" << parsed.total_banks << '\n'
              << "  CondDB ABST0      : " << abst0_val
              << ", UFIELD Bz(0)=" << b0.bz_tesla() << " T\n";
  }

  std::vector<std::string> event_jsons;
  std::vector<delphipp::pipeline::EventRecoOutput> reco_outputs;
  event_jsons.reserve(parsed.events.size());
  reco_outputs.reserve(parsed.events.size());

  for (std::size_t ev_idx = 0; ev_idx < parsed.events.size(); ++ev_idx) {
    auto res = pipeline.reconstruct_event(
        parsed.events[ev_idx],
        opt.run_number,
        static_cast<std::int32_t>(ev_idx + 1));

    if (csv_ofs.is_open()) {
      delphipp::dst::write_csv_event(csv_ofs, res.dst_event);
    }
    if (hm_ofs.is_open()) {
      delphipp::dst::write_hepmc3_event(hm_ofs, res.dst_event);
    }

    event_jsons.push_back(res.tanagra.to_bank_diff_json(res.run, res.event));

    if (opt.verbose) {
      std::cout << "  Event #" << (ev_idx + 1)
                << " (run=" << res.run << ", evt=" << res.event
                << "): banks=" << res.zebra_banks
                << ", raw_blocklets=" << res.raw_blocklets
                << ", hof_hits=" << res.hof_hits
                << ", tracks=" << res.tracks_after_ambiguity
                << "/" << res.tracks_before_ambiguity
                << ", dst_particles=" << res.dst_event.particles.size()
                << ", E_vis=" << res.dst_event.visible_energy_gev << " GeV\n";
    }
    reco_outputs.push_back(std::move(res));
  }

  if (hm_ofs.is_open()) {
    delphipp::dst::write_hepmc3_footer(hm_ofs);
  }

  if (!opt.summary_json_path.empty()) {
    std::ofstream sofs(opt.summary_json_path);
    sofs << "{\n"
         << "  \"input\": \"" << opt.input_fz_path << "\",\n"
         << "  \"detected_zebra_fz\": " << (parsed.detected_zebra_fz ? "true" : "false") << ",\n"
         << "  \"endianness\": \"" << (parsed.is_big_endian ? "big" : "little") << "\",\n"
         << "  \"physical_record_words\": " << parsed.physical_record_words << ",\n"
         << "  \"events_count\": " << reco_outputs.size() << ",\n"
         << "  \"events\": [\n";
    for (std::size_t i = 0; i < reco_outputs.size(); ++i) {
      const auto& r = reco_outputs[i];
      std::string trimmed_top = r.top_bank;
      while (!trimmed_top.empty() && trimmed_top.back() == ' ') {
        trimmed_top.pop_back();
      }
      const float vis_e = std::isfinite(r.dst_event.visible_energy_gev)
                              ? r.dst_event.visible_energy_gev
                              : 0.0f;
      sofs << "    {\"index\": " << (i + 1)
           << ", \"run\": " << r.run
           << ", \"event\": " << r.event
           << ", \"top_bank\": \"" << trimmed_top << "\""
           << ", \"zebra_banks\": " << r.zebra_banks
           << ", \"raw_blocklets\": " << r.raw_blocklets
           << ", \"hof_hits\": " << r.hof_hits
           << ", \"tdr\": " << r.tanagra.tdr_count()
           << ", \"ter\": " << r.tanagra.ter_count()
           << ", \"tkr\": " << r.tanagra.tkr_count()
           << ", \"dst_particles\": " << r.dst_event.particles.size()
           << ", \"visible_energy_gev\": " << vis_e
           << ", \"detectors\": {";
      bool first_det = true;
      for (const auto& [dname, cnt] : r.detector_blocklet_counts) {
        if (!first_det) sofs << ", ";
        first_det = false;
        sofs << "\"" << dname << "\": " << cnt;
      }
      sofs << "}}";
      if (i + 1 < reco_outputs.size()) sofs << ',';
      sofs << '\n';
    }
    sofs << "  ]\n}\n";
  }

  std::string out_payload;
  if (event_jsons.size() == 1) {
    out_payload = event_jsons.front();
  } else {
    out_payload = "[\n";
    for (std::size_t i = 0; i < event_jsons.size(); ++i) {
      out_payload += "  " + event_jsons[i];
      if (i + 1 < event_jsons.size()) {
        out_payload += ',';
      }
      out_payload += '\n';
    }
    out_payload += "]";
  }

  if (!opt.out_json_path.empty()) {
    std::ofstream ofs(opt.out_json_path);
    ofs << out_payload << '\n';
  } else if (!opt.verbose) {
    std::cout << out_payload << '\n';
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  const CliOptions opt = parse_cli(argc, argv);

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

  delphipp::MagneticField bfield;
  bfield.load_delphi_map();

  if (!opt.input_fz_path.empty()) {
    return run_fz_file_reco(opt, abst0.value_or(0.0f), bfield);
  }
  return run_synthetic_event(opt, abst0.value_or(0.0f), bfield);
}
