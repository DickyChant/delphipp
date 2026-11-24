// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <iostream>
#include <string>

#include "delphipp/tanagra_store.hpp"

int main() {
  using namespace delphipp;

  TanagraStore store;

  const TanagraId td_hof = store.add_tdr(TdrBank{
      .detector_id = 29,
      .submod_id = 12,
      .descriptor = 376,
      .num_points = 1,
      .words_per_point = 8,
      .data = {14.0f, 100.0f, -100.0f, 450.0f},
  });
  const TanagraId td_mub = store.add_tdr(TdrBank{
      .detector_id = 14,
      .submod_id = 1,
      .descriptor = 121,
      .num_points = 1,
      .words_per_point = 7,
      .data = {1.5f, 182.5f},
  });

  const TanagraId te_mub = store.add_ter(TerBank{
      .detector_id = 14,
      .submod_id = 1,
      .descriptor = 125,
      .fwords = {445.0f, 120.0f, 50.0f},
      .td_refs = {td_mub},
  });

  const TanagraId tk = store.add_tkr(TkrBank{
      .charge = -1,
      .perigee_params = {0.01f, 0.2f, 1.57f, 0.5f, 0.02f},
      .te_refs = {te_mub},
  });

  assert(store.tdr_count() == 2);
  assert(store.ter_count() == 1);
  assert(store.tkr_count() == 1);
  assert(store.list_tdr(29).size() == 1 && store.list_tdr(29)[0] == td_hof);
  assert(store.list_tdr(14).size() == 1 && store.list_tdr(14)[0] == td_mub);
  assert(store.list_ter(14).size() == 1 && store.list_ter(14)[0] == te_mub);
  assert(store.get_tkr(tk) != nullptr && store.get_tkr(tk)->charge == -1);
  assert(store.get_tdr(999) == nullptr);

  const std::string json = store.to_bank_diff_json(47014, 1);
  assert(json.find("\"run\":47014") != std::string::npos);
  assert(json.find("\"idh\":\"TDR\"") != std::string::npos);
  assert(json.find("\"idh\":\"TER\"") != std::string::npos);

  store.clear();
  assert(store.tdr_count() == 0 && store.ter_count() == 0 && store.tkr_count() == 0);

  std::cout << "test_tanagra_store: PASS\n";
  return 0;
}
