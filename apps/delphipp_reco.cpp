// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// DELPHI++ pure-C++20 reconstruction driver.
// Supports both:
//   1. Self-contained synthetic Z -> mu+ mu- event mode (when --input is omitted).
//   2. Real CERNLIB ZEBRA FZ Exchange file mode (--input <file.fadsim|.sdst|.fz>):
//      parses physical/logical ZEBRA bank structures, populates ZebraStoreView,
//      runs raw detector unpacking (HOF/EGRA, raw blocklets) & Short-DST PA/PV
//      reconstruction, and exports bank_diff JSON, CSV, and HepMC3 Asciiv3.

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
#include "delphipp/io/raw_blocklet.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/reco/hof.hpp"
#include "delphipp/reco/mub.hpp"
#include "delphipp/reco/tof.hpp"
#include "delphipp/tanagra_store.hpp"
#include "delphipp/zebra_store.hpp"

namespace {

struct CliOptions {
  std::string input_fz_path;
  std::string out_json_path;
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

delphipp::PerigeeTrack track_from_p3_charge(float px,
                                            float py,
                                            float pz,
                                            std::int32_t charge,
                                            float bz,
                                            float eps = 0.0f,
                                            float z0 = 0.0f) {
  const float p = std::sqrt(px * px + py * py + pz * pz);
  const float theta = (p > 1e-6f) ? std::acos(std::clamp(pz / p, -1.0f, 1.0f)) : 1.5707963f;
  const float phi0 = std::atan2(py, px);
  const float signed_q = (charge == 0) ? 1.0f : static_cast<float>(charge);
  const float inv_r = (p > 1e-5f) ? (signed_q * bz / p) : 0.0f;
  return delphipp::PerigeeTrack{
      .eps = eps,
      .z0 = z0,
      .theta = theta,
      .phi0 = phi0,
      .inv_r = inv_r,
  };
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

  // Also synthesize corresponding DstEvent if CSV or HepMC3 outputs requested
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
  event_jsons.reserve(parsed.events.size());
  delphipp::reco::hof::HofCalibration hof_calib{};
  delphipp::reco::hof::HofStats hof_stats{};
  hof_stats.reset_run();

  for (std::size_t ev_idx = 0; ev_idx < parsed.events.size(); ++ev_idx) {
    const auto& ev = parsed.events[ev_idx];
    const std::int32_t run = (ev.run_number > 0) ? ev.run_number : opt.run_number;
    const std::int32_t evt = (ev.event_number > 0)
                                 ? ev.event_number
                                 : static_cast<std::int32_t>(ev_idx + 1);

    delphipp::OwnedZebraStore zstore(65536);
    auto zview = zstore.view();
    const auto bank_heads = ev.populate_zebra_store(zview);

    delphipp::TanagraStore tstore;

    // 1. Run HOF raw blocklet unpacking if EGRA/EQRA banks are present
    std::int32_t hof_hits_count = 0;
    if (auto it = bank_heads.find("EGRA"); it != bank_heads.end()) {
      hof_stats.begin_event();
      const std::int32_t mraw =
          bank_heads.contains("EQRA") ? bank_heads.at("EQRA") : 0;
      const auto hof_res = delphipp::reco::hof::unpack_td(
          zview, it->second, mraw, hof_calib, hof_stats);
      hof_hits_count = hof_res.counter_hits;
      if (!hof_res.hits.empty()) {
        std::vector<float> flat;
        flat.reserve(hof_res.hits.size() * 8);
        for (const auto& h : hof_res.hits) {
          flat.insert(flat.end(), h.words.begin(), h.words.end());
        }
        tstore.add_tdr(delphipp::TdrBank{
            .detector_id = delphipp::reco::hof::kDetectorId,
            .submod_id = 1,
            .descriptor = 290,
            .num_points = static_cast<std::int32_t>(hof_res.hits.size()),
            .words_per_point = 8,
            .data = std::move(flat),
        });
      }
    }

    // 2. Summarize raw blocklets from any present DELPHI *RA banks into TanagraStore
    std::size_t raw_blocklets_decoded = 0;
    for (const std::string raw_id : {"EGRA", "EQRA", "EDRA", "EMRA", "EFRA", "RAW "}) {
      if (auto it = bank_heads.find(raw_id); it != bank_heads.end()) {
        const auto blocklets = delphipp::io::read_raw_blocklet_chain(zview, it->second);
        raw_blocklets_decoded += blocklets.size();
        for (const auto& bl : blocklets) {
          tstore.add_tdr(delphipp::TdrBank{
              .detector_id = bl.detector_id,
              .submod_id = bl.subdet_id,
              .descriptor = bl.flags,
              .num_points = 1,
              .words_per_point = static_cast<std::int32_t>(
                  std::min<std::size_t>(bl.data.size(), 4)),
              .data = {
                  static_cast<float>(bl.num_words),
                  static_cast<float>(bl.detector_id),
                  static_cast<float>(bl.subdet_id),
                  abst0_val,
              },
          });
        }
      }
    }

    // 3. Extract PV vertex and PA charged/neutral particles (from .sdst/.fadana)
    delphipp::reco::tracking::BeamspotConstraint beamspot{};
    if (const auto* pv = ev.find_first("PV  "); pv != nullptr && pv->fwords.size() >= 7) {
      const float vx = pv->fwords[4];
      const float vy = pv->fwords[5];
      const float vz = pv->fwords[6];
      if (std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) &&
          std::fabs(vx) < 10.0f && std::fabs(vy) < 10.0f && std::fabs(vz) < 50.0f) {
        beamspot.xyz_cm = {vx, vy, vz};
      }
    }

    std::vector<delphipp::dst::TrackPidInput> charged_tracks;
    std::vector<delphipp::dst::NeutralInput> neutrals;
    const auto pa_banks = ev.find_all("PA  ");
    for (const auto* pa : pa_banks) {
      if (pa->iwords.size() < 20) {
        continue;
      }
      // pxds33.car: IP1 = LPA + 5 + IQ(LPA+4) -> 0-based index 4 + iwords[3]
      const std::int32_t nmod = pa->iwords[3];
      if (nmod < 1 || nmod > 32) {
        continue;
      }
      const std::size_t ip1 = static_cast<std::size_t>(4 + nmod);
      if (ip1 + 10 >= pa->fwords.size()) {
        continue;
      }
      const float px = pa->fwords[ip1 + 2];
      const float py = pa->fwords[ip1 + 3];
      const float pz = pa->fwords[ip1 + 4];
      const float energy = pa->fwords[ip1 + 5];
      const float ch_code_f = pa->fwords[ip1 + 7];
      const float rp = (ip1 + 10 < pa->fwords.size()) ? pa->fwords[ip1 + 10] : 0.0f;
      const float dz = (ip1 + 19 < pa->fwords.size()) ? pa->fwords[ip1 + 19] : 0.0f;

      if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
          !std::isfinite(energy) || energy < 0.0f || energy > 200.0f) {
        continue;
      }
      const int ch_code = static_cast<int>(std::lround(ch_code_f));
      if (ch_code == 1 || ch_code == 2) {
        const std::int32_t q = (ch_code == 1) ? +1 : -1;
        const auto trk = track_from_p3_charge(
            px, py, pz, q, b0.bz,
            std::isfinite(rp) ? rp : 0.0f,
            std::isfinite(dz) ? dz : 0.0f);
        charged_tracks.push_back(delphipp::dst::TrackPidInput{
            .track = trk,
            .charge = q,
        });
        tstore.add_tkr(delphipp::TkrBank{
            .charge = q,
            .perigee_params = {trk.eps, trk.z0, trk.theta, trk.phi0, trk.inv_r},
        });
      } else {
        const float p_mag = std::sqrt(px * px + py * py + pz * pz);
        if (p_mag > 1e-4f) {
          const float scale = 200.0f / p_mag;
          neutrals.push_back(delphipp::dst::NeutralInput{
              .xyz_cm = {px * scale, py * scale, pz * scale},
              .energy_gev = energy,
              .is_ecal = true,
          });
        }
      }
    }

    const auto dst_ev = delphipp::dst::build_short_dst(
        run, evt, b0.bz, charged_tracks, neutrals, beamspot);

    if (csv_ofs.is_open()) {
      delphipp::dst::write_csv_event(csv_ofs, dst_ev);
    }
    if (hm_ofs.is_open()) {
      delphipp::dst::write_hepmc3_event(hm_ofs, dst_ev);
    }

    event_jsons.push_back(tstore.to_bank_diff_json(run, evt));

    if (opt.verbose) {
      std::cout << "  Event #" << (ev_idx + 1) << " (run=" << run << ", evt=" << evt
                << "): banks=" << ev.banks.size()
                << ", raw_blocklets=" << raw_blocklets_decoded
                << ", hof_hits=" << hof_hits_count
                << ", PA_charged=" << charged_tracks.size()
                << ", PA_neutral=" << neutrals.size()
                << ", dst_particles=" << dst_ev.particles.size()
                << ", E_vis=" << dst_ev.visible_energy_gev << " GeV\n";
    }
  }

  if (hm_ofs.is_open()) {
    delphipp::dst::write_hepmc3_footer(hm_ofs);
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

  // Condition DB (94D TPC ABST0 correction)
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
