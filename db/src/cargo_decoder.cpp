// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/db/cargo_decoder.hpp"

#include <algorithm>
#include <cctype>

namespace delphipp::db {

namespace {

// KAPACK character code table (0..99, printable ASCII subset starting at ' ')
int char_to_ka_code(char c) noexcept {
  const auto u = static_cast<unsigned char>(c);
  if (u >= 32 && u <= 126) {
    return static_cast<int>(u - 31);  // 1..95
  }
  return 0;
}

char ka_code_to_char(int code) noexcept {
  if (code >= 1 && code <= 95) {
    return static_cast<char>(code + 31);
  }
  return ' ';
}

}  // namespace

std::vector<std::int32_t> pack_kapack_key(std::string_view key, int max_pack) {
  const int mp = std::max(1, std::min(max_pack, 4));
  std::vector<std::int32_t> words;
  std::size_t i = 0;
  while (i < key.size()) {
    std::int32_t acc = 0;
    for (int k = 0; k < mp; ++k) {
      acc *= 100;
      if (i < key.size()) {
        acc += char_to_ka_code(key[i++]);
      }
    }
    words.push_back(acc);
  }
  return words;
}

std::string unpack_kapack_key(std::span<const std::int32_t> words, int max_pack) {
  const int mp = std::max(1, std::min(max_pack, 4));
  std::string out;
  out.reserve(words.size() * static_cast<std::size_t>(mp));
  for (std::int32_t w : words) {
    std::string chunk(static_cast<std::size_t>(mp), ' ');
    std::int32_t v = (w < 0) ? -w : w;
    for (int k = mp - 1; k >= 0; --k) {
      const int code = v % 100;
      v /= 100;
      chunk[static_cast<std::size_t>(k)] = (code == 0) ? '\0' : ka_code_to_char(code);
    }
    for (const char c : chunk) {
      if (c != '\0') out.push_back(c);
    }
  }
  while (!out.empty() && out.back() == ' ') {
    out.pop_back();
  }
  return out;
}

}  // namespace delphipp::db
