// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <sstream>

#include "delphipp/db/cargo_decoder.hpp"
#include "delphipp/dst/dst_writer.hpp"
#include "delphipp/dst/short_dst.hpp"
#include "delphipp/magnetic_field.hpp"

int main() {
  using namespace delphipp;

  // 1. Real DELPHI measured UFIELD map (41 barrel + 39 endcap layers from ufield.car)
  MagneticField field;
  field.load_delphi_map();
  const auto b_center = field.evaluate({0.0f, 0.0f, 0.0f});
  assert(std::fabs(b_center.bz_tesla() - 1.2f) < 1e-3f);

  // Inside HAC barrel return yoke (r = 350 cm, z = 100 cm): return flux reverses Bz!
  const auto b_yoke = field.evaluate({350.0f, 0.0f, 100.0f});
  assert(b_yoke.bmag > 0.0f);
  assert(b_yoke.bz_tesla() < 0.0f);  // Iron return yoke opposes central solenoid

  // Endcap return yoke (r = 150 cm, z = 400 cm)
  const auto b_endcap = field.evaluate({150.0f, 0.0f, 400.0f});
  assert(b_endcap.bmag > 0.0f);

  // SAT superconducting quadrupole (CSAT = 0.0007609)
  const auto b_sat = field.evaluate({3.0f, 4.0f, 450.0f});
  assert(std::fabs(b_sat.bx - 0.0007609f * 4.0f) < 1e-6f);
  assert(std::fabs(b_sat.by - 0.0007609f * 3.0f) < 1e-6f);

  // 2. CARGO base-100 key pack/unpack round-trip
  const std::string key = "GLOB$TPC*.ENP0";
  const auto packed = db::pack_kapack_key(key, 4);
  const auto unpacked = db::unpack_kapack_key(packed, 4);
  assert(unpacked == key);

  // 3. Columnar CSV, NDJSON, and HepMC3-ASCII DST writers
  std::vector<dst::TrackPidInput> trks = {
      {.track = {.eps = 0.0f, .z0 = 0.0f, .theta = std::numbers::pi_v<float> * 0.5f,
                 .phi0 = 0.5f, .inv_r = b_center.bz / 45.6f},
       .charge = -1, .has_muon_hit = true},
  };
  const auto ev = dst::build_short_dst(47014, 1, b_center.bz, trks, {});

  std::ostringstream csv;
  dst::write_csv_header(csv);
  dst::write_csv_event(csv, ev);
  assert(csv.str().find("47014,1,1,-1,13,") != std::string::npos);

  std::ostringstream hepmc;
  dst::write_hepmc3_header(hepmc);
  dst::write_hepmc3_event(hepmc, ev);
  dst::write_hepmc3_footer(hepmc);
  assert(hepmc.str().find("HepMC::Version 3.02.06") != std::string::npos);
  assert(hepmc.str().find("P 1 -1 13 ") != std::string::npos);

  std::cout << "test_ufield_cargo_dst: PASS\n";
  return 0;
}
