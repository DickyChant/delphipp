// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/io/raw_blocklet.hpp"

namespace delphipp::io {

std::vector<RawBlocklet> read_raw_blocklet_chain(const ZebraStoreView& store,
                                                 std::int32_t first_link) {
  std::vector<RawBlocklet> out;
  std::int32_t loc = first_link;
  std::size_t guard = 0;
  constexpr std::size_t kMaxChain = 10'000;

  while (loc > 0 && guard++ < kMaxChain) {
    const std::int32_t nd = store.bank_nd(loc);
    if (nd < 2) {
      break;
    }
    const auto hdr2 = static_cast<std::uint32_t>(store.iq(loc + 2));
    RawBlocklet blk;
    blk.bank_link = loc;
    blk.detector_id = static_cast<std::int32_t>((hdr2 >> 8) & 0xFFu);
    blk.subdet_id = static_cast<std::int32_t>((hdr2 >> 16) & 0x0Fu);
    blk.flags = static_cast<std::int32_t>((hdr2 >> 20) & 0xFFFu);
    blk.num_words = nd;
    blk.data.reserve(static_cast<std::size_t>(nd));
    for (std::int32_t k = 1; k <= nd; ++k) {
      blk.data.push_back(store.iq(loc + k));
    }
    out.push_back(std::move(blk));
    loc = store.bank_next(loc);
  }
  return out;
}

}  // namespace delphipp::io
