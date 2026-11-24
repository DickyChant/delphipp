// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/zebra_store.hpp"

#include <cstring>

namespace delphipp {

ZebraStoreView::ZebraStoreView(std::span<std::int32_t> lq_words)
    : lq_(lq_words) {
  if (lq_.size() <= fortran::kLqToIqOffset) {
    throw std::invalid_argument("ZebraStoreView requires more than 8 words");
  }
}

void ZebraStoreView::check_lq_index(std::int32_t idx) const {
  if (idx < 1 || static_cast<std::size_t>(idx) > lq_.size()) {
    throw std::out_of_range("ZEBRA LQ index out of bounds: " +
                            std::to_string(idx));
  }
}

std::int32_t ZebraStoreView::lq(std::int32_t idx) const {
  check_lq_index(idx);
  return lq_[static_cast<std::size_t>(idx - 1)];
}

void ZebraStoreView::set_lq(std::int32_t idx, std::int32_t value) {
  check_lq_index(idx);
  lq_[static_cast<std::size_t>(idx - 1)] = value;
}

std::int32_t ZebraStoreView::iq(std::int32_t idx) const {
  return lq(idx + static_cast<std::int32_t>(fortran::kLqToIqOffset));
}

void ZebraStoreView::set_iq(std::int32_t idx, std::int32_t value) {
  set_lq(idx + static_cast<std::int32_t>(fortran::kLqToIqOffset), value);
}

float ZebraStoreView::q(std::int32_t idx) const {
  return std::bit_cast<float>(iq(idx));
}

void ZebraStoreView::set_q(std::int32_t idx, float value) {
  set_iq(idx, std::bit_cast<std::int32_t>(value));
}

std::string ZebraStoreView::bank_idh(std::int32_t l) const {
  const std::int32_t raw = iq(l - 4);
  char chars[4];
  std::memcpy(chars, &raw, 4);
  return std::string(chars, 4);
}

OwnedZebraStore::OwnedZebraStore(std::size_t n_words)
    : storage_(n_words + fortran::kLqToIqOffset, 0) {}

}  // namespace delphipp
