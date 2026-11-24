// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "delphipp/reco/ambiguity_solver.hpp"

int main() {
  using namespace delphipp::reco::tracking;

  // Track 1 (score 10) shares TE 101 with Track 2 (score 6) and TE 102 with Track 3 (score 6).
  // Optimal conflict-free solution is {Track 2, Track 3} (score 12) + independent Track 4 (score 3) = 15.
  std::vector<CandidateTrack> cands = {
      {.id = 1, .te_ids = {101, 102}, .score = 10.0},
      {.id = 2, .te_ids = {101, 201}, .score = 6.0},
      {.id = 3, .te_ids = {102, 301}, .score = 6.0},
      {.id = 4, .te_ids = {401}, .score = 3.0},
  };

  const auto res = solve_ambiguities(cands);
  assert(!res.hit_node_limit);
  assert(std::fabs(res.total_score - 15.0) < 1e-9);
  assert(res.accepted_track_ids.size() == 3);

  auto ids = res.accepted_track_ids;
  std::sort(ids.begin(), ids.end());
  assert((ids == std::vector<std::int32_t>{2, 3, 4}));

  // Deterministic node budget cutoff
  const auto limited = solve_ambiguities(cands, SolverConfig{.max_nodes = 2});
  assert(limited.hit_node_limit);
  assert(limited.nodes_visited == 2);

  std::cout << "test_ambiguity_solver: PASS\n";
  return 0;
}
