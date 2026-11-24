// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 in-memory TANAGRA bank store (TDR, TER, TKR) replacing CERNLIB/ZEBRA
// TANAGRA calls (TDETEC, TELEM, TLIST, TGET) for DELPHI++ reconstruction modules.

#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace delphipp {

using TanagraId = std::int32_t;

// Detector Data Bank (TD/TDR created via TDETEC)
struct TdrBank {
  TanagraId id = 0;
  std::int32_t detector_id = 0;
  std::int32_t submod_id = 0;
  std::int32_t descriptor = 0;
  std::int32_t num_points = 0;
  std::int32_t words_per_point = 0;
  std::array<float, 3> header_floats{};
  std::vector<float> data;
};

// Track Element Bank (TE/TER created via TELEM)
struct TerBank {
  TanagraId id = 0;
  std::int32_t detector_id = 0;
  std::int32_t submod_id = 0;
  std::int32_t descriptor = 0;
  std::vector<std::int32_t> iwords;  // 1-based words 1..7 typically integer
  std::vector<float> fwords;         // 1-based words 8..N typically float
  std::vector<TanagraId> td_refs;    // referenced TDR bank IDs
};

// Reconstructed Track Bank (TK/TKR)
struct TkrBank {
  TanagraId id = 0;
  std::int32_t charge = 0;
  std::array<float, 5> perigee_params{};  // (eps, z, theta, phi, 1/R)
  std::vector<TanagraId> te_refs;
};

class TanagraStore {
 public:
  void clear() noexcept;

  // Mirrors TDETEC: registers a TDR bank and returns its assigned 1-based ID.
  TanagraId add_tdr(TdrBank bank);

  // Mirrors TELEM: registers a TER bank and returns its assigned 1-based ID.
  TanagraId add_ter(TerBank bank);

  // Registers a TKR bank and returns its assigned 1-based ID.
  TanagraId add_tkr(TkrBank bank);

  // Mirrors TLIST('A', 'TDR'/'TER', 'DET', det_id, ...)
  [[nodiscard]] std::vector<TanagraId> list_tdr(std::int32_t detector_id) const;
  [[nodiscard]] std::vector<TanagraId> list_ter(std::int32_t detector_id) const;

  // Mirrors TGET('TDR'/'TER', id, ...)
  [[nodiscard]] const TdrBank* get_tdr(TanagraId id) const noexcept;
  [[nodiscard]] const TerBank* get_ter(TanagraId id) const noexcept;
  [[nodiscard]] const TkrBank* get_tkr(TanagraId id) const noexcept;

  [[nodiscard]] std::size_t tdr_count() const noexcept { return tdrs_.size(); }
  [[nodiscard]] std::size_t ter_count() const noexcept { return ters_.size(); }
  [[nodiscard]] std::size_t tkr_count() const noexcept { return tkrs_.size(); }

  // Export event banks in the JSON schema expected by tools/bank_diff.py.
  [[nodiscard]] std::string to_bank_diff_json(std::int32_t run,
                                              std::int32_t event) const;

 private:
  std::vector<TdrBank> tdrs_;
  std::vector<TerBank> ters_;
  std::vector<TkrBank> tkrs_;
};

}  // namespace delphipp
