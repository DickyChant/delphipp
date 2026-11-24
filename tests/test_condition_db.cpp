// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

#include "delphipp/db/condition_db.hpp"

int main() {
  using namespace delphipp::db;

  ConditionDatabase db;

  // Base full record before 1994 (drift velocity word 1 = 6.72, ABST0 word 14 = 1.147)
  db.add_record(ConditionRecord{
      .table = "CALB",
      .path = "GLOB$TPC*.GLOB$ENP0",
      .is_word_update = false,
      .valid_from = {931101, 0},
      .entered_on = {940601, 120000},
      .field = "SLOW",
      .words = {{1, 6.72f}, {14, 1.147f}},
  });

  // Exact format of simana_v94d/ddb/tpct0_94c.asc
  std::istringstream asc(R"(
*DBRECORD
CALB   GLOB$TPC*.GLOB$ENP0.W
  940101,000001,990101,120000
*DBFIELD SLOW
!94C reprocessing: absolute T0 1.147 - 0.04492 time bins (TPC z +-232 um)
 3,SLOW,14., 1.10208
*DBRECORD
CALB   GLOB$TPC*.GLOB$ENP0.W
  950101,000001,990101,120000
*DBFIELD SLOW
!94C reprocessing: restore absolute T0 for 1995
 3,SLOW,14., 1.147
*DBRECORD
CALB   GLOB$TPC*.GLOB$ENP0.W
  940101,000001,260101,120000
*DBFIELD SLOW
!Post-2000 entry should be excluded by DELANA default cutoff (2000-01-01)
 3,SLOW,14., 9.999
)");

  const std::size_t parsed = db.load_ascii(asc);
  assert(parsed == 3);

  // Before 1994: ABST0 = 1.147
  const auto pre94 = db.resolve_field("CALB", "GLOB$TPC*.GLOB$ENP0", "SLOW", {931115, 120000});
  assert(std::fabs(pre94.at(1) - 6.72f) < 1e-5f);
  assert(std::fabs(pre94.at(14) - 1.147f) < 1e-5f);

  // During 1994 (run 47014, May 1994): drift velocity untouched, ABST0 = 1.10208
  const auto in94 = db.resolve_field("CALB", "GLOB$TPC*.GLOB$ENP0", "SLOW", {940526, 120000});
  assert(std::fabs(in94.at(1) - 6.72f) < 1e-5f);
  assert(std::fabs(in94.at(14) - 1.10208f) < 1e-5f);

  // In 1995: ABST0 restored to 1.147
  const auto in95 = db.resolve_word("CALB", "GLOB$TPC*.GLOB$ENP0", "SLOW", 14, {950201, 0});
  assert(in95.has_value());
  assert(std::fabs(*in95 - 1.147f) < 1e-5f);

  std::cout << "test_condition_db: PASS\n";
  return 0;
}
