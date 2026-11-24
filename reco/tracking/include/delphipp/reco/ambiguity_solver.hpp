// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 global track ambiguity solver inspired by DELPHI's fxsolv.c / fxambi.car.
// Finds a maximum-score conflict-free set of candidate tracks (no shared TE IDs)
// using deterministic branch-and-bound with a fixed node-visit budget.

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/tanagra_store.hpp"

namespace delphipp::reco::tracking {

struct CandidateTrack {
  std::int32_t id = 0;
  std::vector<TanagraId> te_ids;  // Track Elements used by this candidate
  double score = 0.0;             // Higher is better (e.g. detector weight + chi2 prob)
};

struct SolverConfig {
  // Deterministic node budget (replaces wall-clock timeout in legacy fxsolv.c)
  std::uint64_t max_nodes = 100'000;
};

struct SolverResult {
  std::vector<std::int32_t> accepted_track_ids;
  double total_score = 0.0;
  std::uint64_t nodes_visited = 0;
  bool hit_node_limit = false;
};

// Select the highest-scoring conflict-free subset of `candidates`
// (two candidates conflict if they share any TE ID).
[[nodiscard]] SolverResult solve_ambiguities(
    std::span<const CandidateTrack> candidates,
    const SolverConfig& cfg = {});

}  // namespace delphipp::reco::tracking
