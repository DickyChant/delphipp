// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/io/fz_scanner.hpp"

#include <array>
#include <cstring>
#include <unordered_set>

#include "delphipp/io/fz_exchange.hpp"

namespace delphipp::io {

namespace {

std::uint32_t read_u32_be(const unsigned char* p) noexcept {
  return (static_cast<std::uint32_t>(p[0]) << 24) |
         (static_cast<std::uint32_t>(p[1]) << 16) |
         (static_cast<std::uint32_t>(p[2]) << 8) |
         (static_cast<std::uint32_t>(p[3]));
}

std::uint32_t read_u32_le(const unsigned char* p) noexcept {
  return (static_cast<std::uint32_t>(p[0])) |
         (static_cast<std::uint32_t>(p[1]) << 8) |
         (static_cast<std::uint32_t>(p[2]) << 16) |
         (static_cast<std::uint32_t>(p[3]) << 24);
}

bool match_sync(const unsigned char* p, bool big_endian) noexcept {
  auto rd = big_endian ? read_u32_be : read_u32_le;
  for (std::size_t i = 0; i < 4; ++i) {
    if (rd(p + i * 4) != kFzSyncWords[i]) {
      return false;
    }
  }
  return true;
}

const std::unordered_set<std::string>& known_delphi_bank_ids() {
  static const std::unordered_set<std::string> kIds = {
      "HEAD", "DFIX", "LDTO", "PV  ", "PA  ", "TK  ", "TS  ", "TE  ", "TD  ",
      "EGRA", "EQRA", "EDRA", "EMRA", "EFRA", "MAIN", "EMNC", "MUID", "ELID",
      "HAID", "BTAG", "SIMU", "EVHD", "RUNH",
  };
  return kIds;
}

}  // namespace

FzScanSummary scan_fz_stream(std::istream& in) {
  FzScanSummary sum;
  std::vector<unsigned char> buf(
      (std::istreambuf_iterator<char>(in)),
      std::istreambuf_iterator<char>());
  sum.file_bytes = buf.size();
  if (buf.size() < 32) {
    return sum;
  }

  // Detect ZEBRA FZ sync words at start (or first steering block) in BE or LE
  bool big_endian = true;
  std::size_t first_off = buf.size();
  for (std::size_t off = 0; off + 32 <= buf.size(); off += 4) {
    if (match_sync(buf.data() + off, true)) {
      big_endian = true;
      first_off = off;
      break;
    }
    if (match_sync(buf.data() + off, false)) {
      big_endian = false;
      first_off = off;
      break;
    }
  }

  if (first_off >= buf.size()) {
    return sum;
  }

  sum.detected_zebra_fz = true;
  sum.is_big_endian = big_endian;
  auto rd = big_endian ? read_u32_be : read_u32_le;

  // Word 4 of steering block holds physical record length in 32-bit words (typically 900)
  const std::uint32_t raw_len = rd(buf.data() + first_off + 16);
  const std::uint32_t phys_words = (raw_len & 0x00FFFFFFu);
  sum.physical_record_words = (phys_words >= 16 && phys_words <= 65536) ? phys_words : 900u;
  const std::size_t phys_bytes = static_cast<std::size_t>(sum.physical_record_words) * 4u;

  for (std::size_t off = first_off; off + 32 <= buf.size(); off += phys_bytes) {
    if (match_sync(buf.data() + off, big_endian)) {
      ++sum.steering_blocks;
      // Word 7 (offset 28) is pilot offset within physical block (in words)
      const std::uint32_t pilot_off = rd(buf.data() + off + 28) & 0xFFFFu;
      if (pilot_off >= 8 && (off + pilot_off * 4 + 8) <= buf.size()) {
        const std::uint32_t lr_len = rd(buf.data() + off + pilot_off * 4);
        const std::uint32_t lr_type = rd(buf.data() + off + pilot_off * 4 + 4);
        if (lr_len > 0 && lr_type <= 8) {
          ++sum.logical_records;
          ++sum.record_type_counts[lr_type];
        }
      }
    } else {
      ++sum.fast_blocks;
    }
  }

  // Scan aligned 4-byte words for known DELPHI Hollerith bank IDs
  const auto& bank_set = known_delphi_bank_ids();
  for (std::size_t off = 0; off + 4 <= buf.size(); off += 4) {
    std::string id(reinterpret_cast<const char*>(buf.data() + off), 4);
    if (bank_set.contains(id)) {
      ++sum.bank_id_counts[id];
    }
  }

  return sum;
}

}  // namespace delphipp::io
