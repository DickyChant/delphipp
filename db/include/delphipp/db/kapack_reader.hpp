// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 reader for CERNLIB KAPACK v2 binary direct-access files
// (DELPHI condition database snapshots FOR050–FOR056, see kapack600.car
// KAXMAK & KAXGET).

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <optional>
#include <string>
#include <vector>

namespace delphipp::db {

inline constexpr std::size_t kKapackRootControlWords = 24;  // NBCW1 in kapack600.car

// Decoded KAPACK Block 1 header (matches IA(1..24) in kapack600.car KAXMAK)
struct KapackFileHeader {
  std::int32_t block_number = 0;   // IA(1)  == 1
  std::int32_t free_blocks = 0;    // IA(2)
  std::int32_t root_block = 0;     // IA(4)
  std::int32_t tree_levels = 0;    // IA(5)
  std::int32_t first_entry = 0;    // IA(6)  == NBCW1 + 1 (25)
  std::int32_t next_free_word = 0; // IA(7)
  std::int32_t version = 0;        // IA(13) == 2 (KAVRSN)
  std::int32_t num_blocks = 0;     // IA(14) == NB
  std::int32_t block_words = 0;    // IA(15) == LBLK (e.g. 2048 words = 8192 bytes)
  std::int32_t max_pack = 0;       // IA(16) == MAXPAK (e.g. 4)
  std::uint32_t file_id = 0;       // IA(24) == KAXFID()
};

struct KapackBlock {
  std::int32_t block_number = 0;
  std::int32_t tree_level = 0;     // IA(5): 1 = leaf data block, >1 = index block
  std::int32_t first_entry = 0;    // IA(6)
  std::int32_t next_free_word = 0; // IA(7)
  std::vector<std::int32_t> words;
};

struct KapackEntry {
  std::vector<std::int32_t> key_words;      // Packed key words (excluding length word)
  std::vector<std::int32_t> payload_words;  // Child block pointer (level>1) or record data (level==1)
};

class KapackReader {
 public:
  // Opens and validates a KAPACK binary stream (e.g. decompressed FOR050..FOR056).
  explicit KapackReader(std::istream& in);

  [[nodiscard]] const KapackFileHeader& header() const noexcept { return header_; }
  [[nodiscard]] std::size_t block_bytes() const noexcept {
    return static_cast<std::size_t>(header_.block_words) * 4u;
  }

  // Reads 1-based direct-access block `block_num` (1 <= block_num <= num_blocks).
  [[nodiscard]] std::optional<KapackBlock> read_block(std::int32_t block_num);

  // Parses all entries in a single KAPACK block (following KAXGET IA(6)..IA(7) chain).
  [[nodiscard]] static std::vector<KapackEntry> parse_block_entries(
      const KapackBlock& blk);

  // Traverses the KAPACK B-tree from `header().root_block` down to level-1 leaf blocks
  // and returns all leaf entries (up to `max_entries`).
  [[nodiscard]] std::vector<KapackEntry> list_leaf_entries(
      std::size_t max_entries = 100'000);

 private:
  std::istream& in_;
  KapackFileHeader header_{};
};

}  // namespace delphipp::db
