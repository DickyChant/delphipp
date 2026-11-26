// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#include "delphipp/io/fz_zebra_parser.hpp"
#include "delphipp/magnetic_field.hpp"
#include "delphipp/pipeline/event_pipeline.hpp"

namespace {

std::int32_t ascii_to_word(const char s[4]) {
  return (static_cast<std::int32_t>(static_cast<unsigned char>(s[0])) << 24) |
         (static_cast<std::int32_t>(static_cast<unsigned char>(s[1])) << 16) |
         (static_cast<std::int32_t>(static_cast<unsigned char>(s[2])) << 8) |
         (static_cast<std::int32_t>(static_cast<unsigned char>(s[3])));
}

delphipp::io::ParsedZebraBank make_pa(std::int32_t tk_id,
                                      float px,
                                      float py,
                                      float pz,
                                      float energy,
                                      float charge_code) {
  delphipp::io::ParsedZebraBank b;
  b.idh = "PA  ";
  b.idh_raw = ascii_to_word("PA  ");
  b.nd = 32;
  b.iwords.assign(32, 0);
  b.fwords.assign(32, 0.0f);
  b.iwords[0] = tk_id;
  b.iwords[3] = 1;
  b.iwords[4] = 26;
  auto set_f = [&](std::size_t idx, float v) {
    b.fwords[idx] = v;
    b.iwords[idx] = std::bit_cast<std::int32_t>(v);
  };
  set_f(5 + 2, px);
  set_f(5 + 3, py);
  set_f(5 + 4, pz);
  set_f(5 + 5, energy);
  set_f(5 + 7, charge_code);
  return b;
}

delphipp::io::ParsedZebraBank make_bank(const char idh[4],
                                        std::int32_t orig_addr,
                                        std::int32_t nl,
                                        std::int32_t ns,
                                        const std::vector<std::int32_t>& links,
                                        const std::vector<std::int32_t>& words) {
  delphipp::io::ParsedZebraBank b;
  b.idh = idh;
  b.idh_raw = ascii_to_word(idh);
  b.orig_addr = orig_addr;
  b.nl = nl;
  b.ns = ns;
  b.nd = static_cast<std::int32_t>(words.size());
  b.links = links;
  b.iwords = words;
  for (auto w : words) b.fwords.push_back(std::bit_cast<float>(w));
  return b;
}

}  // namespace

int main() {
  using namespace delphipp;

  MagneticField bfield;
  bfield.load_delphi_map();
  pipeline::EventReconstructionPipeline pipe(bfield, 1.10208f);

  // 1. Relocated DELSIM tree across all 18 subdetectors: TOP -> RAW -> PMBs -> RAWD
  {
    io::ParsedFzEvent raw_ev;
    raw_ev.run_number = -47014;
    raw_ev.event_number = 7;
    raw_ev.has_header_event_number = true;
    raw_ev.top_bank = "TOP ";

    // Tracking hits aligned around phi ~ 0.52..0.65 rad + calorimeters + muon + RICH
    auto b_vd  = make_bank("RAWD", 100, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((1u << 23) | (2u << 18) | (100u << 8) | 40u),
         static_cast<std::int32_t>((1u << 25) | (1u << 23) | (2u << 18) | (100u << 8) | 40u)});
    auto p_vd  = make_bank("RVD ", 110, 5, 5, {0, 0, 0, 0, 100}, {15, 0});

    auto b_id  = make_bank("RAWD", 120, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((2u << 15) | (2u << 10) | 100u),
         static_cast<std::int32_t>((2u << 15) | (10u << 10) | 120u)});
    auto p_id  = make_bank("RID ", 130, 5, 5, {0, 0, 0, 0, 120}, {15, 0});

    auto b_tpc = make_bank("RAWD", 140, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((0u << 24) | (0u << 20) | (25u << 8) | 45u),
         static_cast<std::int32_t>((0u << 24) | (8u << 20) | (20u << 8) | 50u)});
    auto p_tpc = make_bank("RTPC", 150, 5, 5, {0, 0, 0, 0, 140}, {15, 0});

    auto b_od  = make_bank("RAWD", 160, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((2u << 19) | (0u << 16) | (20u << 8) | 20u),
         static_cast<std::int32_t>((2u << 19) | (2u << 16) | (22u << 8) | 20u)});
    auto p_od  = make_bank("ROD ", 170, 5, 5, {0, 0, 0, 0, 160}, {15, 0});

    auto b_fca = make_bank("RAWD", 180, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((0u << 10) | 600u),
         static_cast<std::int32_t>((1u << 10) | 600u)});
    auto p_fca = make_bank("RFCA", 190, 5, 5, {0, 0, 0, 0, 180}, {15, 0});

    auto b_fcb = make_bank("RAWD", 200, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((0u << 10) | 620u),
         static_cast<std::int32_t>((1u << 10) | 620u)});
    auto p_fcb = make_bank("RFCB", 210, 5, 5, {0, 0, 0, 0, 200}, {15, 0});

    auto b_hpc = make_bank("RAWD", 220, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((2u << 21) | (2u << 17) | (40u << 8) | 80u), 0});
    auto p_hpc = make_bank("RHPC", 230, 5, 5, {0, 0, 0, 0, 220}, {15, 0});

    auto b_emf = make_bank("RAWD", 240, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((40u << 18) | (40u << 12) | 400u), 0});
    auto p_emf = make_bank("REMF", 250, 5, 5, {0, 0, 0, 0, 240}, {15, 0});

    // RHAC far in phi (port=10 -> neutral hadron in DstEvent)
    auto b_hac = make_bank("RAWD", 260, 0, 0, {},
        {7, 0, 0x5A5A, 0, (10 << 16) | (10 << 8) | 20, (10 << 16) | (11 << 8) | 10, 0});
    auto p_hac = make_bank("RHAC", 270, 5, 5, {0, 0, 0, 0, 260}, {15, 0});

    auto b_sat = make_bank("RAWD", 280, 0, 0, {},
        {5, 0, 0, 0, static_cast<std::int32_t>((12u << 15) | (1u << 12) | 40u)});
    auto p_sat = make_bank("RSAT", 290, 5, 5, {0, 0, 0, 0, 280}, {15, 0});

    auto b_lum = make_bank("RAWD", 300, 0, 0, {},
        {5, 0, 0, 0, static_cast<std::int32_t>((13u << 15) | (1u << 12) | 40u)});
    auto p_lum = make_bank("RLUM", 310, 5, 5, {0, 0, 0, 0, 300}, {15, 0});

    auto b_tof = make_bank("RAWD", 320, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((3u << 27) | (4u << 17) | 200u), 0});
    auto p_tof = make_bank("RTOF", 330, 5, 5, {0, 0, 0, 0, 320}, {15, 0});

    auto b_hof = make_bank("RAWD", 340, 0, 0, {},
        {8, 0, (22 << 8) | (1 << 16), 0, 0, 1, 0, 0, 0});
    auto p_hof = make_bank("RHOF", 350, 5, 5, {0, 0, 0, 0, 340}, {15, 0});

    auto b_mub = make_bank("RAWD", 360, 0, 0, {},
        {5, 0, 0, 0, static_cast<std::int32_t>((2u << 22) | 100u)});
    auto p_mub = make_bank("RMUB", 370, 5, 5, {0, 0, 0, 0, 360}, {15, 0});

    auto b_muf = make_bank("RAWD", 380, 0, 0, {},
        {5, 0, 0, 0, static_cast<std::int32_t>((0u << 10) | 50u)});
    auto p_muf = make_bank("RMUF", 390, 5, 5, {0, 0, 0, 0, 380}, {15, 0});

    auto b_rib = make_bank("RAWD", 400, 0, 0, {},
        {5, 0, 0, 0, static_cast<std::int32_t>((2u << 10) | 200u)});
    auto p_rib = make_bank("RRIB", 410, 5, 5, {0, 0, 0, 0, 400}, {15, 0});

    auto b_rif = make_bank("RAWD", 420, 0, 0, {},
        {5, 0, 0, 0, static_cast<std::int32_t>((2u << 10) | 200u)});
    auto p_rif = make_bank("RRIF", 430, 5, 5, {0, 0, 0, 0, 420}, {15, 0});

    // Legacy EDRA for MUS surround muon chambers
    auto b_edra = make_bank("EDRA", 440, 0, 0, {},
        {6, 0, 0, 0, static_cast<std::int32_t>((0u << 14) | (5u << 8) | 60u),
         static_cast<std::int32_t>((1u << 14) | (6u << 8) | 60u)});

    // RAW (links 6..21 -> PMBs)
    std::vector<std::int32_t> raw_links(22, 0);
    raw_links[6 - 1]  = 110;  // RVD
    raw_links[7 - 1]  = 130;  // RID
    raw_links[8 - 1]  = 150;  // RTPC
    raw_links[9 - 1]  = 410;  // RRIB
    raw_links[10 - 1] = 170;  // ROD
    raw_links[11 - 1] = 230;  // RHPC
    raw_links[12 - 1] = 330;  // RTOF
    raw_links[13 - 1] = 270;  // RHAC
    raw_links[14 - 1] = 370;  // RMUB
    raw_links[15 - 1] = 190;  // RFCA
    raw_links[16 - 1] = 430;  // RRIF
    raw_links[17 - 1] = 290;  // RSAT
    raw_links[18 - 1] = 210;  // RFCB
    raw_links[19 - 1] = 250;  // REMF
    raw_links[20 - 1] = 350;  // RHOF
    raw_links[21 - 1] = 390;  // RMUF
    raw_links[22 - 1] = 310;  // RLUM
    auto raw_bank = make_bank("RAW ", 480, 22, 22, raw_links, {15, 0});
    auto top_bank = make_bank("TOP ", 500, 7, 7, {480, 0, 0, 0, 0, 0, 0}, {1});

    raw_ev.banks = {
        b_vd, p_vd, b_id, p_id, b_tpc, p_tpc, b_od, p_od, b_fca, p_fca,
        b_fcb, p_fcb, b_hpc, p_hpc, b_emf, p_emf, b_hac, p_hac, b_sat, p_sat,
        b_lum, p_lum, b_tof, p_tof, b_hof, p_hof, b_mub, p_mub, b_muf, p_muf,
        b_rib, p_rib, b_rif, p_rif, b_edra, raw_bank, top_bank,
    };

    const auto out = pipe.reconstruct_event(raw_ev);
    assert(out.run == -47014);
    assert(out.event == 7);
    assert(out.hof_hits == 1);
    assert(out.raw_blocklets == 18);
    for (const char* det : {"VD", "ID", "TPC", "OD", "FCA", "FCB", "HPC", "EMF",
                            "HAC", "SAT", "LUM", "TOF", "HOF", "MUB", "MUF",
                            "MUS", "RIB", "RIF"}) {
      assert(out.detector_blocklet_counts.at(det) == 1);
    }
    assert(out.kalman_refits >= 1);
    assert(out.tracks_after_ambiguity >= 1);
    assert(out.tanagra.tkr_count() >= 1);
    bool has_charged = false;
    bool has_neutral = false;
    for (const auto& ptcl : out.dst_event.particles) {
      if (ptcl.charge != 0) has_charged = true;
      if (ptcl.charge == 0) has_neutral = true;
    }
    assert(has_charged && has_neutral);
    assert(std::isfinite(out.dst_event.visible_energy_gev));
    assert(out.dst_event.visible_energy_gev > 5.0f);
  }

  // 2. DST event: ambiguity conflict between two charged PA candidates sharing TK=77 + 1 neutral
  {
    io::ParsedFzEvent dst_ev;
    dst_ev.run_number = 47014;
    dst_ev.event_number = 8;

    dst_ev.banks.push_back(make_pa(77, 10.0f, 5.0f, 10.0f, 15.0f, 1.0f));
    dst_ev.banks.push_back(make_pa(77, 1.0f, 2.0f, 2.0f, 3.0f, 1.0f));
    dst_ev.banks.push_back(make_pa(88, -2.0f, 0.0f, 4.0f, 4.47f, 0.0f));

    const auto out = pipe.reconstruct_event(dst_ev);
    assert(out.tracks_before_ambiguity == 2);
    assert(out.tracks_after_ambiguity == 1);
    assert(out.tanagra.tkr_count() == 1);
    assert(out.dst_event.particles.size() == 2);
    assert(out.dst_event.visible_energy_gev > 19.0f);
  }

  // 3. DELSIM SIM_ truth record (SP vertex + ST tracks)
  {
    io::ParsedFzEvent sim_ev;
    sim_ev.run_number = -47014;
    sim_ev.event_number = 3;
    sim_ev.has_header_event_number = true;
    sim_ev.top_bank = "SIM_";

    std::vector<std::int32_t> sp_w(6, 0);
    auto set_fw = [](std::vector<std::int32_t>& w, std::size_t idx, float v) {
      w[idx] = std::bit_cast<std::int32_t>(v);
    };
    set_fw(sp_w, 1, 0.01f);
    set_fw(sp_w, 2, -0.02f);
    set_fw(sp_w, 3, 0.35f);
    sim_ev.banks.push_back(make_bank("SP  ", 0, 2, 1, {}, sp_w));

    std::vector<std::int32_t> st1(13, 0);
    set_fw(st1, 1, 20.0f);
    set_fw(st1, 2, 5.0f);
    set_fw(st1, 3, 12.0f);
    set_fw(st1, 4, 23.85f);
    set_fw(st1, 6, -1.0f);
    set_fw(st1, 7, 0.13957f);
    sim_ev.banks.push_back(make_bank("ST  ", 0, 5, 1, {}, st1));

    std::vector<std::int32_t> st_beam(13, 0);
    set_fw(st_beam, 1, 0.0f);
    set_fw(st_beam, 2, 0.0f);
    set_fw(st_beam, 3, 45.625f);
    set_fw(st_beam, 4, 45.625f);
    set_fw(st_beam, 6, -1.0f);
    set_fw(st_beam, 7, 0.000511f);
    sim_ev.banks.push_back(make_bank("ST  ", 0, 5, 1, {}, st_beam));

    const auto out = pipe.reconstruct_event(sim_ev);
    assert(out.dst_event.particles.size() == 1);
    assert(out.dst_event.particles[0].charge == -1);
    assert(std::isfinite(out.dst_event.visible_energy_gev));
    assert(std::fabs(out.dst_event.visible_energy_gev - 23.85f) < 0.5f);
    const std::string dj = out.dst_event.to_json();
    assert(dj.find("inf") == std::string::npos);
    assert(dj.find("nan") == std::string::npos);
  }

  std::cout << "test_reco_pipeline: PASS\n";
  return 0;
}
