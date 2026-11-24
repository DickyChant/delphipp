// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 reader for DELPHI standard FASTBUS/LTD raw-data ZEBRA blocklets
// (used across HOF, MUS, MUB, MUF, TOF, TPC, HPC raw banks in ZebraStoreView).

#pragma once

#include <cstdint>
#include <vector>

#include "delphipp/zebra_store.hpp"

namespace delphipp::io {

struct RawBlocklet {
  std::int32_t bank_link = 0;     // 1-based ZEBRA link L
  std::int32_t detector_id = 0;   // Bits 8..15 of IQ(L+2)
  std::int32_t subdet_id = 0;     // Bits 16..19 of IQ(L+2)
  std::int32_t flags = 0;         // Bits 20..31 of IQ(L+2)
  std::int32_t num_words = 0;     // IQ(L-1)
  std::vector<std::int32_t> data; // IQ(L+1 .. L+nd)
};

// Walk a linear ZEBRA chain starting at `first_link` and decode each raw blocklet.
[[nodiscard]] std::vector<RawBlocklet> read_raw_blocklet_chain(
    const ZebraStoreView& store,
    std::int32_t first_link);

}  // namespace delphipp::io
