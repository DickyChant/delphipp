// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/db/kapack_reader.hpp"

#include <array>
#include <stdexcept>
#include <unordered_set>

namespace delphipp::db {

namespace {

std::int32_t read_le32(const unsigned char* p) noexcept {
  const std::uint32_t u =
      static_cast<std::uint32_t>(p[0]) |
      (static_cast<std::uint32_t>(p[1]) << 8) |
      (static_cast<std::uint32_t>(p[2]) << 16) |
      (static_cast<std::uint32_t>(p[3]) << 24);
  return static_cast<std::int32_t>(u);
}

}  // namespace

KapackReader::KapackReader(std::istream& in) : in_(in) {
  in_.seekg(0, std::ios::beg);
  std::array<unsigned char, kKapackRootControlWords * 4> raw{};
  if (!in_.read(reinterpret_cast<char*>(raw.data()),
                static_cast<std::streamsize>(raw.size()))) {
    throw std::runtime_error("KAPACK: truncated root block header");
  }

  auto w = [&](std::size_t one_based) -> std::int32_t {
    return read_le32(raw.data() + (one_based - 1) * 4);
  };

  header_.block_number = w(1);
  header_.free_blocks = w(2);
  header_.root_block = w(4);
  header_.tree_levels = w(5);
  header_.first_entry = w(6);
  header_.next_free_word = w(7);
  header_.version = w(13);
  header_.num_blocks = w(14);
  header_.block_words = w(15);
  header_.max_pack = w(16);
  header_.file_id = static_cast<std::uint32_t>(w(24));

  if (header_.block_number != 1 || header_.version != 2 ||
      header_.block_words < static_cast<std::int32_t>(kKapackRootControlWords) ||
      header_.num_blocks < 1) {
    throw std::runtime_error("KAPACK: invalid root block header");
  }
}

std::optional<KapackBlock> KapackReader::read_block(std::int32_t block_num) {
  if (block_num < 1 || block_num > header_.num_blocks) {
    return std::nullopt;
  }
  const std::size_t nbytes = block_bytes();
  const std::streamoff offset =
      static_cast<std::streamoff>(static_cast<std::size_t>(block_num - 1) * nbytes);
  in_.clear();
  in_.seekg(offset, std::ios::beg);

  std::vector<unsigned char> raw(nbytes);
  if (!in_.read(reinterpret_cast<char*>(raw.data()),
                static_cast<std::streamsize>(nbytes))) {
    return std::nullopt;
  }

  KapackBlock blk;
  blk.words.resize(static_cast<std::size_t>(header_.block_words));
  for (std::size_t i = 0; i < blk.words.size(); ++i) {
    blk.words[i] = read_le32(raw.data() + i * 4);
  }
  blk.block_number = blk.words[0];
  blk.tree_level = blk.words[4];
  blk.first_entry = blk.words[5];
  blk.next_free_word = blk.words[6];
  return blk;
}

std::vector<KapackEntry> KapackReader::parse_block_entries(const KapackBlock& blk) {
  std::vector<KapackEntry> out;
  const auto nwords = static_cast<std::int32_t>(blk.words.size());
  std::int32_t loc = blk.first_entry;  // 1-based (IA(6))
  const std::int32_t limit = std::min(blk.next_free_word, nwords + 1);

  if (loc <= 6 || loc >= limit) {
    return out;
  }

  // Matches KAXGET: IA(LOC) = entry length, IA(LOC+1) = key length (incl. length word)
  while (loc < limit) {
    const std::int32_t entry_len = blk.words[static_cast<std::size_t>(loc - 1)];
    if (entry_len < 2 || loc + entry_len > limit) {
      break;
    }
    const std::int32_t key_len = blk.words[static_cast<std::size_t>(loc)];
    if (key_len < 1 || 1 + key_len > entry_len) {
      break;
    }

    KapackEntry entry;
    for (std::int32_t k = 1; k < key_len; ++k) {
      entry.key_words.push_back(blk.words[static_cast<std::size_t>(loc + k)]);
    }
    for (std::int32_t p = 1 + key_len; p < entry_len; ++p) {
      entry.payload_words.push_back(blk.words[static_cast<std::size_t>(loc - 1 + p)]);
    }
    out.push_back(std::move(entry));
    loc += entry_len;
  }
  return out;
}

std::vector<KapackEntry> KapackReader::list_leaf_entries(std::size_t max_entries) {
  std::vector<KapackEntry> leaves;
  std::unordered_set<std::int32_t> visited;

  auto walk = [&](auto&& self, std::int32_t block_num) -> void {
    if (leaves.size() >= max_entries || !visited.insert(block_num).second) {
      return;
    }
    const auto blk = read_block(block_num);
    if (!blk.has_value()) {
      return;
    }
    auto entries = parse_block_entries(*blk);
    if (blk->tree_level <= 1) {
      for (auto& e : entries) {
        if (leaves.size() >= max_entries) break;
        leaves.push_back(std::move(e));
      }
    } else {
      for (const auto& e : entries) {
        if (!e.payload_words.empty()) {
          self(self, e.payload_words.front());
        }
      }
    }
  };

  walk(walk, header_.root_block);
  return leaves;
}

}  // namespace delphipp::db
