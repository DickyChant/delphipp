// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 reader for CERNLIB KAPACK v2 binary direct-access files
// (DELPHI condition database snapshots FOR050–FOR056, see kapack600.car KAXMAK).

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
  std::vector<std::int32_t> words;
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

 private:
  std::istream& in_;
  KapackFileHeader header_{};
};

}  // namespace delphipp::db
