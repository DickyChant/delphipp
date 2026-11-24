// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// C++20 view over a ZEBRA store (either owned by C++ or backed by Fortran
// COMMON /UXCOM/ during the strangler-pattern transition).

#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "delphipp/ux_common.hpp"

namespace delphipp {

// Standard ZEBRA bank header offsets relative to bank address L (1-based in Fortran):
//   LQ(L)     : link to next bank in linear chain (next)
//   LQ(L + 1) : link to up bank (up)
//   LQ(L + 2) : link to origin link (origin)
//   IQ(L - 5) : numeric bank ID (idn)
//   IQ(L - 4) : Hollerith 4-char bank ID (idh)
//   IQ(L - 3) : total number of links (nl)
//   IQ(L - 2) : number of structural links (ns)
//   IQ(L - 1) : number of data words (nd)
//   IQ(L)     : status word (status)
//   IQ(L + k) / Q(L + k) : data word k (1 <= k <= nd)
//   LQ(L - j) : structural/reference link j (1 <= j <= nl)
//
// Because IQ(1) == LQ(9), IQ(idx) is at lq_base[idx + 7] (0-based),
// and LQ(idx) is at lq_base[idx - 1] (0-based).

class ZebraStoreView {
 public:
  // Construct a non-owning view over a ZEBRA store starting at LQ(1).
  // `words` must contain at least `kLqToIqOffset` words.
  explicit ZebraStoreView(std::span<std::int32_t> lq_words);

  [[nodiscard]] std::size_t size_words() const noexcept { return lq_.size(); }

  // 1-based Fortran LQ(idx) access.
  [[nodiscard]] std::int32_t lq(std::int32_t idx) const;
  void set_lq(std::int32_t idx, std::int32_t value);

  // 1-based Fortran IQ(idx) access (equivalent to LQ(idx + 8)).
  [[nodiscard]] std::int32_t iq(std::int32_t idx) const;
  void set_iq(std::int32_t idx, std::int32_t value);

  // 1-based Fortran Q(idx) access (bit-cast float view of IQ(idx)).
  [[nodiscard]] float q(std::int32_t idx) const;
  void set_q(std::int32_t idx, float value);

  // Read standard ZEBRA bank metadata for bank at 1-based link `l`.
  [[nodiscard]] std::int32_t bank_next(std::int32_t l) const { return lq(l); }
  [[nodiscard]] std::int32_t bank_up(std::int32_t l) const { return lq(l + 1); }
  [[nodiscard]] std::int32_t bank_idn(std::int32_t l) const { return iq(l - 5); }
  [[nodiscard]] std::string bank_idh(std::int32_t l) const;
  [[nodiscard]] std::int32_t bank_nl(std::int32_t l) const { return iq(l - 3); }
  [[nodiscard]] std::int32_t bank_ns(std::int32_t l) const { return iq(l - 2); }
  [[nodiscard]] std::int32_t bank_nd(std::int32_t l) const { return iq(l - 1); }
  [[nodiscard]] std::int32_t bank_status(std::int32_t l) const { return iq(l); }

  // Structural/reference link j (1 <= j <= nl) of bank at `l`.
  [[nodiscard]] std::int32_t bank_link(std::int32_t l, std::int32_t j) const {
    return lq(l - j);
  }

  // Data word k (1 <= k <= nd) of bank at `l`.
  [[nodiscard]] std::int32_t bank_idata(std::int32_t l, std::int32_t k) const {
    return iq(l + k);
  }
  [[nodiscard]] float bank_fdata(std::int32_t l, std::int32_t k) const {
    return q(l + k);
  }

 private:
  std::span<std::int32_t> lq_;

  void check_lq_index(std::int32_t idx) const;
};

// Owning ZEBRA store for standalone C++ tests and future Fortran-free execution.
class OwnedZebraStore {
 public:
  explicit OwnedZebraStore(std::size_t n_words = fortran::kDefaultNSizeQ);

  [[nodiscard]] ZebraStoreView view() noexcept { return ZebraStoreView(storage_); }
  [[nodiscard]] std::span<std::int32_t> raw_words() noexcept { return storage_; }

 private:
  std::vector<std::int32_t> storage_;
};

}  // namespace delphipp
