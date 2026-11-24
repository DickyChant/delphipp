// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/tanagra_store.hpp"

#include <iomanip>
#include <sstream>

namespace delphipp {

void TanagraStore::clear() noexcept {
  tdrs_.clear();
  ters_.clear();
  tkrs_.clear();
}

TanagraId TanagraStore::add_tdr(TdrBank bank) {
  bank.id = static_cast<TanagraId>(tdrs_.size() + 1);
  tdrs_.push_back(std::move(bank));
  return tdrs_.back().id;
}

TanagraId TanagraStore::add_ter(TerBank bank) {
  bank.id = static_cast<TanagraId>(ters_.size() + 1);
  ters_.push_back(std::move(bank));
  return ters_.back().id;
}

TanagraId TanagraStore::add_tkr(TkrBank bank) {
  bank.id = static_cast<TanagraId>(tkrs_.size() + 1);
  tkrs_.push_back(std::move(bank));
  return tkrs_.back().id;
}

std::vector<TanagraId> TanagraStore::list_tdr(std::int32_t detector_id) const {
  std::vector<TanagraId> out;
  for (const auto& b : tdrs_) {
    if (b.detector_id == detector_id) {
      out.push_back(b.id);
    }
  }
  return out;
}

std::vector<TanagraId> TanagraStore::list_ter(std::int32_t detector_id) const {
  std::vector<TanagraId> out;
  for (const auto& b : ters_) {
    if (b.detector_id == detector_id) {
      out.push_back(b.id);
    }
  }
  return out;
}

const TdrBank* TanagraStore::get_tdr(TanagraId id) const noexcept {
  if (id < 1 || static_cast<std::size_t>(id) > tdrs_.size()) {
    return nullptr;
  }
  return &tdrs_[static_cast<std::size_t>(id - 1)];
}

const TerBank* TanagraStore::get_ter(TanagraId id) const noexcept {
  if (id < 1 || static_cast<std::size_t>(id) > ters_.size()) {
    return nullptr;
  }
  return &ters_[static_cast<std::size_t>(id - 1)];
}

const TkrBank* TanagraStore::get_tkr(TanagraId id) const noexcept {
  if (id < 1 || static_cast<std::size_t>(id) > tkrs_.size()) {
    return nullptr;
  }
  return &tkrs_[static_cast<std::size_t>(id - 1)];
}

std::string TanagraStore::to_bank_diff_json(std::int32_t run,
                                            std::int32_t event) const {
  std::ostringstream os;
  os << std::setprecision(7);
  os << "{\"events\":[{\"run\":" << run << ",\"event\":" << event << ",\"banks\":[";
  bool first_bank = true;

  auto comma = [&]() {
    if (!first_bank) {
      os << ',';
    }
    first_bank = false;
  };

  for (const auto& b : tdrs_) {
    comma();
    os << "{\"idh\":\"TDR\",\"idn\":" << b.id
       << ",\"iwords\":{\"1\":" << b.detector_id
       << ",\"2\":" << b.submod_id
       << ",\"3\":" << b.descriptor
       << ",\"4\":" << b.num_points
       << ",\"5\":" << b.words_per_point << "},\"fwords\":{";
    for (std::size_t i = 0; i < b.data.size(); ++i) {
      if (i > 0) os << ',';
      os << '"' << (i + 1) << "\":" << b.data[i];
    }
    os << "}}";
  }

  for (const auto& b : ters_) {
    comma();
    os << "{\"idh\":\"TER\",\"idn\":" << b.id
       << ",\"iwords\":{\"1\":" << b.detector_id
       << ",\"2\":" << b.submod_id
       << ",\"4\":" << b.descriptor << "},\"fwords\":{";
    for (std::size_t i = 0; i < b.fwords.size(); ++i) {
      if (i > 0) os << ',';
      os << '"' << (i + 1) << "\":" << b.fwords[i];
    }
    os << "}}";
  }

  os << "]}]}";
  return os.str();
}

}  // namespace delphipp
