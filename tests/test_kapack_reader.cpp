// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <vector>

#include "delphipp/db/kapack_reader.hpp"

namespace {

void write_le32(std::vector<unsigned char>& buf, std::size_t one_based_word,
                std::uint32_t v) {
  const std::size_t off = (one_based_word - 1) * 4;
  buf[off + 0] = static_cast<unsigned char>(v & 0xFFu);
  buf[off + 1] = static_cast<unsigned char>((v >> 8) & 0xFFu);
  buf[off + 2] = static_cast<unsigned char>((v >> 16) & 0xFFu);
  buf[off + 3] = static_cast<unsigned char>((v >> 24) & 0xFFu);
}

}  // namespace

int main() {
  using namespace delphipp::db;

  // Build a 3-block KAPACK v2 stream matching FOR050's control-word layout
  constexpr std::size_t kBlockWords = 64;
  constexpr std::size_t kNumBlocks = 3;
  std::vector<unsigned char> img(kNumBlocks * kBlockWords * 4, 0);

  // Block 1 root header (matches IA(1..24) in kapack600.car KAXMAK)
  write_le32(img, 1, 1u);             // block_number
  write_le32(img, 2, 1u);             // free_blocks
  write_le32(img, 4, 1u);             // root_block
  write_le32(img, 5, 1u);             // tree_levels
  write_le32(img, 6, 25u);            // first_entry = NBCW1 + 1
  write_le32(img, 7, 47u);            // next_free_word
  write_le32(img, 13, 2u);            // KAVRSN = 2
  write_le32(img, 14, kNumBlocks);    // NB
  write_le32(img, 15, kBlockWords);   // LBLK
  write_le32(img, 16, 4u);            // MAXPAK
  write_le32(img, 24, 0x01E648F1u);   // file_id from FOR050

  // Block 2 sentinel
  write_le32(img, kBlockWords + 1, 2u);
  write_le32(img, kBlockWords + 10, 47014u);

  std::string data(reinterpret_cast<const char*>(img.data()), img.size());
  std::istringstream iss(data, std::ios::in | std::ios::binary);

  KapackReader reader(iss);
  assert(reader.header().version == 2);
  assert(reader.header().num_blocks == static_cast<std::int32_t>(kNumBlocks));
  assert(reader.header().block_words == static_cast<std::int32_t>(kBlockWords));
  assert(reader.header().max_pack == 4);
  assert(reader.header().file_id == 0x01E648F1u);

  const auto blk2 = reader.read_block(2);
  assert(blk2.has_value());
  assert(blk2->block_number == 2);
  assert(blk2->words[9] == 47014);
  assert(!reader.read_block(99).has_value());

  std::cout << "test_kapack_reader: PASS\n";
  return 0;
}
