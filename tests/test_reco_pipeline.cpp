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

}  // namespace

int main() {
  using namespace delphipp;

  MagneticField bfield;
  bfield.load_delphi_map();
  pipeline::EventReconstructionPipeline pipe(bfield, 1.10208f);

  // 1. Raw multi-detector event: EGRA (HOF simulated hit) + EMRA (MUB) + TPRA (TPC)
  {
    io::ParsedFzEvent raw_ev;
    raw_ev.run_number = 47014;
    raw_ev.event_number = 7;

    io::ParsedZebraBank egra;
    egra.idh = "EGRA";
    egra.idh_raw = ascii_to_word("EGRA");
    egra.iwords = {0, (22 << 8) | (1 << 16), 0, 0, 1, 0, 0, 0};
    egra.nd = static_cast<std::int32_t>(egra.iwords.size());
    for (auto w : egra.iwords) egra.fwords.push_back(std::bit_cast<float>(w));

    io::ParsedZebraBank emra;
    emra.idh = "EMRA";
    emra.idh_raw = ascii_to_word("EMRA");
    emra.iwords = {0, (14 << 8) | (2 << 16), 42, 99};
    emra.nd = static_cast<std::int32_t>(emra.iwords.size());
    for (auto w : emra.iwords) emra.fwords.push_back(std::bit_cast<float>(w));

    io::ParsedZebraBank tpra;
    tpra.idh = "TPRA";
    tpra.idh_raw = ascii_to_word("TPRA");
    tpra.iwords = {0, (5 << 8) | (1 << 16), 123};
    tpra.nd = static_cast<std::int32_t>(tpra.iwords.size());
    for (auto w : tpra.iwords) tpra.fwords.push_back(std::bit_cast<float>(w));

    raw_ev.banks = {egra, emra, tpra};

    const auto out = pipe.reconstruct_event(raw_ev);
    assert(out.run == 47014);
    assert(out.event == 7);
    assert(out.hof_hits == 1);
    assert(out.raw_blocklets == 3);
    assert(out.detector_blocklet_counts.at("HOF") == 1);
    assert(out.detector_blocklet_counts.at("MUB") == 1);
    assert(out.detector_blocklet_counts.at("TPC") == 1);
    assert(out.tanagra.tdr_count() >= 4);
    assert(out.tanagra.ter_count() == 3);
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

  std::cout << "test_reco_pipeline: PASS\n";
  return 0;
}
