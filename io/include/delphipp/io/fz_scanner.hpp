// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Real CERNLIB ZEBRA FZ Exchange binary scanner.
// Detects ZEBRA steering blocks (sync words 0x0123CDEF, 0x80708070,
// 0x4321ABCD, 0x80618061) in either big-endian or little-endian byte order,
// parses physical record headers & fast blocks, and extracts logical record
// statistics and 4-char Hollerith bank IDs from real DELSIM/DELANA .fz/.sdst/.fadsim files.

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <map>
#include <string>
#include <vector>

namespace delphipp::io {

struct FzScanSummary {
  bool detected_zebra_fz = false;
  bool is_big_endian = true;
  std::size_t file_bytes = 0;
  std::uint32_t physical_record_words = 0;
  std::size_t steering_blocks = 0;
  std::size_t fast_blocks = 0;
  std::size_t logical_records = 0;
  std::map<std::uint32_t, std::size_t> record_type_counts;
  std::map<std::string, std::size_t> bank_id_counts;
};

// Scan a binary stream for real CERNLIB ZEBRA FZ Exchange steering blocks,
// logical records, and known DELPHI Hollerith bank IDs.
[[nodiscard]] FzScanSummary scan_fz_stream(std::istream& in);

}  // namespace delphipp::io
