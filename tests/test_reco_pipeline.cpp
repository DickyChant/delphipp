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

  // 1. Relocated DELSIM tree: TOP -> RAW -> RHAC / RTOF / RTPC / RHOF -> RAWD
  {
    io::ParsedFzEvent raw_ev;
    raw_ev.run_number = -47014;
    raw_ev.event_number = 7;
    raw_ev.has_header_event_number = true;
    raw_ev.top_bank = "TOP ";

    // RHAC RAWD with 0x5A5A marker + 2 hits (energy = adc * 0.17 GeV -> 20*0.17=3.4 GeV)
    auto rawd_hac = make_bank("RAWD", 100, 0, 0, {},
                              {7, 0, 0x5A5A, 0, (1 << 16) | (10 << 8) | 20,
                               (1 << 16) | (11 << 8) | 10, 0});
    auto rhac = make_bank("RHAC", 120, 5, 5, {0, 0, 0, 0, 100}, {15, 0});

    // RTOF RAWD with FASTBUS TDC hit (slot=3, chan=4, tdc=200 -> 100 ns)
    auto rawd_tof = make_bank("RAWD", 140, 0, 0, {},
                              {6, 0, 0, 0, static_cast<std::int32_t>((3u << 27) | (4u << 17) | 200u), 0});
    auto rtof = make_bank("RTOF", 160, 5, 5, {0, 0, 0, 0, 140}, {15, 0});

    // RTPC RAWD with FADC pad hit (endplate=0, sector=2, row=3, time=25, charge=40)
    auto rawd_tpc = make_bank("RAWD", 180, 0, 0, {},
                              {6, 0, 0, 0, static_cast<std::int32_t>((1u << 24) | (2u << 20) | (25u << 8) | 40u), 0});
    auto rtpc = make_bank("RTPC", 200, 5, 5, {0, 0, 0, 0, 180}, {15, 0});

    // RHOF RAWD simulated HOF hit
    auto rawd_hof = make_bank("RAWD", 220, 0, 0, {},
                              {8, 0, (22 << 8) | (1 << 16), 0, 0, 1, 0, 0, 0});
    auto rhof = make_bank("RHOF", 240, 5, 5, {0, 0, 0, 0, 220}, {15, 0});

    // RAW (links: 8->RTPC=200, 12->RTOF=160, 13->RHAC=120, 20->RHOF=240)
    std::vector<std::int32_t> raw_links(22, 0);
    raw_links[8 - 1] = 200;
    raw_links[12 - 1] = 160;
    raw_links[13 - 1] = 120;
    raw_links[20 - 1] = 240;
    auto raw_bank = make_bank("RAW ", 280, 22, 22, raw_links, {15, 0});
    auto top_bank = make_bank("TOP ", 320, 7, 7, {280, 0, 0, 0, 0, 0, 0}, {1});

    raw_ev.banks = {rawd_hac, rhac, rawd_tof, rtof, rawd_tpc, rtpc, rawd_hof, rhof, raw_bank, top_bank};

    const auto out = pipe.reconstruct_event(raw_ev);
    assert(out.run == -47014);
    assert(out.event == 7);
    assert(out.hof_hits == 1);
    assert(out.raw_blocklets == 4);
    assert(out.detector_blocklet_counts.at("HAC") == 1);
    assert(out.detector_blocklet_counts.at("TOF") == 1);
    assert(out.detector_blocklet_counts.at("TPC") == 1);
    assert(out.detector_blocklet_counts.at("HOF") == 1);
    assert(out.tanagra.tdr_count() >= 6);
    assert(out.tanagra.ter_count() >= 6);
    assert(!out.dst_event.particles.empty());  // HCAL neutral cluster propagated into ShortDST
    assert(out.dst_event.visible_energy_gev > 4.0f);
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

    const auto out = pipe.reconstruct_event(sim_ev);
    assert(out.dst_event.particles.size() == 1);
    assert(out.dst_event.particles[0].charge == -1);
    assert(out.dst_event.visible_energy_gev > 20.0f);
  }

  std::cout << "test_reco_pipeline: PASS\n";
  return 0;
}
