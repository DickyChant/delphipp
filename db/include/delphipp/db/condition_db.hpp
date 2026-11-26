// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 DELPHI condition database store and ASCII (.asc) update parser.
// Supports CARGO/DDAP validity timestamps, entry-date precedence, and `.W`
// partial word updates (e.g. simana_v94d/ddb/tpct0_94c.asc).

#pragma once

#include <compare>
#include <cstdint>
#include <istream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace delphipp::db {

// DELPHI condition DB timestamp: YYMMDD and HHMMSS (LEP era: 80..99 -> 1980..1999,
// 00..79 -> 2000..2079, matching CARGO/DELANA convention).
struct DbTimestamp {
  std::int32_t yymmdd = 0;
  std::int32_t hhmmss = 0;

  [[nodiscard]] std::int64_t sortable_key() const noexcept {
    const std::int32_t yy = yymmdd / 10000;
    const std::int32_t mmdd = yymmdd % 10000;
    const std::int32_t year = (yy >= 80) ? (1900 + yy) : (2000 + yy);
    return static_cast<std::int64_t>(year * 10000 + mmdd) * 1'000'000LL +
           static_cast<std::int64_t>(hhmmss);
  }

  auto operator<=>(const DbTimestamp& other) const noexcept {
    return sortable_key() <=> other.sortable_key();
  }
  bool operator==(const DbTimestamp& other) const noexcept = default;
};

struct ConditionRecord {
  std::string table;         // e.g. "CALB", "GEOM", "SENS"
  std::string path;          // e.g. "GLOB$TPC*.GLOB$ENP0"
  bool is_word_update = false;  // true for ".W" partial word updates
  DbTimestamp valid_from{};
  DbTimestamp entered_on{};
  std::string field;         // e.g. "SLOW", "SHAP", "CHAN"
  // 1-based word index -> float value
  std::map<std::int32_t, float> words;
};

class ConditionDatabase {
 public:
  // Add a single record directly.
  void add_record(ConditionRecord rec);

  // Parse DELPHI ASCII update stream (*DBRECORD / *DBFIELD format, e.g. tpct0_94c.asc).
  // Returns the number of records parsed.
  std::size_t load_ascii(std::istream& in);

  // Load records from a binary CERNLIB KAPACK v2 stream (e.g. FOR050–FOR056).
  // Returns the number of leaf entries decoded into ConditionRecords.
  std::size_t load_kapack(std::istream& in, std::size_t max_entries = 100'000);

  [[nodiscard]] std::size_t record_count() const noexcept { return records_.size(); }

  // Resolve merged field words for (table, path, field) at `query_time`.
  // Records with `entered_on >= entry_cutoff` are excluded (DELANA excludes
  // entries from 2000-01-01 00:00:00 onwards by default).
  // Applies matching records in ascending (valid_from, entered_on) order;
  // full records replace the field, `.W` records overlay individual words.
  [[nodiscard]] std::map<std::int32_t, float> resolve_field(
      std::string_view table,
      std::string_view path,
      std::string_view field,
      DbTimestamp query_time,
      DbTimestamp entry_cutoff = DbTimestamp{0, 0}) const;

  [[nodiscard]] std::optional<float> resolve_word(
      std::string_view table,
      std::string_view path,
      std::string_view field,
      std::int32_t one_based_word,
      DbTimestamp query_time,
      DbTimestamp entry_cutoff = DbTimestamp{0, 0}) const;

 private:
  std::vector<ConditionRecord> records_;
};

}  // namespace delphipp::db
