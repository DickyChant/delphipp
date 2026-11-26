// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 CERNLIB ZEBRA FZ Exchange bank-tree parser.
// Unpacks physical steering/fast blocks (both big-endian and little-endian),
// reconstructs logical records & ZEBRA bank structures ([IDN, IDH, NL, NS, ND,
// STATUS, links..., data...]), groups them into events delimited by HEAD banks
// or logical EventNormal records, and populates OwnedZebraStore / ZebraStoreView
// chains for downstream DELPHI++ reconstruction modules.

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <map>
#include <string>
#include <vector>

#include "delphipp/zebra_store.hpp"

namespace delphipp::io {

struct ParsedZebraBank {
  std::int32_t idn = 0;
  std::string idh;        // 4-char Hollerith bank ID (e.g. "HEAD", "EGRA", "PA  ", "PV  ")
  std::int32_t idh_raw = 0;
  std::int32_t nl = 0;
  std::int32_t ns = 0;
  std::int32_t nd = 0;
  std::int32_t status = 0;
  std::vector<std::int32_t> links;   // NL link words preceding the bank header if present
  std::vector<std::int32_t> iwords;  // ND integer/raw data words
  std::vector<float> fwords;         // ND IEEE-754 float bit-cast views of iwords
};

struct ParsedFzEvent {
  std::int32_t run_number = 0;
  std::int32_t event_number = 0;
  bool has_header_event_number = false;
  std::string top_bank;
  std::vector<std::int32_t> user_header;
  std::vector<ParsedZebraBank> banks;

  [[nodiscard]] const ParsedZebraBank* find_first(const std::string& idh) const noexcept;
  [[nodiscard]] std::vector<const ParsedZebraBank*> find_all(const std::string& idh) const;

  // Populate `store` with all banks in this event, linking same-ID banks via
  // LQ(L) linear chains, and return the head 1-based address L per bank ID.
  [[nodiscard]] std::map<std::string, std::int32_t> populate_zebra_store(
      ZebraStoreView& store, std::int32_t start_lq_offset = 100) const;
};

struct FzParseResult {
  bool detected_zebra_fz = false;
  bool is_big_endian = true;
  std::size_t file_bytes = 0;
  std::uint32_t physical_record_words = 0;
  std::size_t steering_blocks = 0;
  std::size_t fast_blocks = 0;
  std::size_t total_banks = 0;
  std::vector<ParsedZebraBank> preamble_banks;
  std::vector<ParsedFzEvent> events;
  std::map<std::string, std::size_t> bank_counts;
};

// Parse an entire CERNLIB ZEBRA FZ Exchange binary stream (`.fz`, `.fadsim`,
// `.fadana`, `.sdst`, or DELPHI++ FzExchangeWriter streams) into events and banks.
// Stops cleanly once `max_events > 0` complete events have been parsed.
[[nodiscard]] FzParseResult parse_fz_stream(std::istream& in,
                                            std::size_t max_events = 0);

}  // namespace delphipp::io
