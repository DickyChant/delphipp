// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 CARGO / KAPACK base-100 key packer/unpacker for DELPHI
// condition database binary snapshots (FOR050–FOR056).

#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace delphipp::db {

// Pack an ASCII key string into KAPACK base-100 32-bit words (MAXPAK=4 chars/word).
[[nodiscard]] std::vector<std::int32_t> pack_kapack_key(std::string_view key,
                                                        int max_pack = 4);

// Unpack KAPACK base-100 32-bit key words back into an ASCII string.
[[nodiscard]] std::string unpack_kapack_key(std::span<const std::int32_t> words,
                                            int max_pack = 4);

}  // namespace delphipp::db
