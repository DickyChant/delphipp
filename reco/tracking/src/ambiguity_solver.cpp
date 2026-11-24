// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/reco/ambiguity_solver.hpp"

#include <algorithm>
#include <numeric>
#include <unordered_set>

namespace delphipp::reco::tracking {

namespace {

bool tracks_conflict(const CandidateTrack& a, const CandidateTrack& b) {
  for (const TanagraId te : a.te_ids) {
    if (std::find(b.te_ids.begin(), b.te_ids.end(), te) != b.te_ids.end()) {
      return true;
    }
  }
  return false;
}

}  // namespace

SolverResult solve_ambiguities(std::span<const CandidateTrack> candidates,
                               const SolverConfig& cfg) {
  const std::size_t n = candidates.size();
  SolverResult result;
  if (n == 0) {
    return result;
  }

  // Deterministic ordering by descending score, tie-breaking by ID
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) {
    if (candidates[i].score != candidates[j].score) {
      return candidates[i].score > candidates[j].score;
    }
    return candidates[i].id < candidates[j].id;
  });

  // Precompute conflict adjacency list in sorted order
  std::vector<std::vector<std::size_t>> conflicts(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      if (tracks_conflict(candidates[order[i]], candidates[order[j]])) {
        conflicts[i].push_back(j);
        conflicts[j].push_back(i);
      }
    }
  }

  // Suffix sum of positive scores for admissible upper-bound pruning
  std::vector<double> suffix_max(n + 1, 0.0);
  for (std::size_t i = n; i > 0; --i) {
    const double s = candidates[order[i - 1]].score;
    suffix_max[i - 1] = suffix_max[i] + (s > 0.0 ? s : 0.0);
  }

  std::vector<std::size_t> current_sel;
  std::vector<std::size_t> best_sel;
  std::vector<int> blocked(n, 0);
  double best_score = 0.0;

  auto dfs = [&](auto&& self, std::size_t idx, double cur_score) -> void {
    if (result.nodes_visited >= cfg.max_nodes) {
      result.hit_node_limit = true;
      return;
    }
    ++result.nodes_visited;

    if (cur_score > best_score) {
      best_score = cur_score;
      best_sel = current_sel;
    }
    if (idx >= n || cur_score + suffix_max[idx] <= best_score) {
      return;
    }

    // Branch 1: include candidate `idx` if not blocked and score > 0
    if (blocked[idx] == 0 && candidates[order[idx]].score > 0.0) {
      current_sel.push_back(idx);
      for (const std::size_t nb : conflicts[idx]) {
        ++blocked[nb];
      }
      self(self, idx + 1, cur_score + candidates[order[idx]].score);
      for (const std::size_t nb : conflicts[idx]) {
        --blocked[nb];
      }
      current_sel.pop_back();
    }

    // Branch 2: skip candidate `idx`
    if (!result.hit_node_limit) {
      self(self, idx + 1, cur_score);
    }
  };

  dfs(dfs, 0, 0.0);

  result.total_score = best_score;
  result.accepted_track_ids.reserve(best_sel.size());
  for (const std::size_t idx : best_sel) {
    result.accepted_track_ids.push_back(candidates[order[idx]].id);
  }
  return result;
}

}  // namespace delphipp::reco::tracking
