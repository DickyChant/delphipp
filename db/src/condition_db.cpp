// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/db/condition_db.hpp"

#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "delphipp/db/cargo_decoder.hpp"
#include "delphipp/db/kapack_reader.hpp"

namespace delphipp::db {

namespace {

std::string trim(std::string_view s) {
  std::size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
    ++start;
  }
  std::size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
    --end;
  }
  return std::string(s.substr(start, end - start));
}

std::vector<std::string> split_csv(std::string_view line) {
  std::vector<std::string> parts;
  std::size_t start = 0;
  while (start <= line.size()) {
    const std::size_t comma = line.find(',', start);
    if (comma == std::string_view::npos) {
      parts.push_back(trim(line.substr(start)));
      break;
    }
    parts.push_back(trim(line.substr(start, comma - start)));
    start = comma + 1;
  }
  return parts;
}

}  // namespace

void ConditionDatabase::add_record(ConditionRecord rec) {
  records_.push_back(std::move(rec));
}

std::size_t ConditionDatabase::load_ascii(std::istream& in) {
  std::size_t added = 0;
  std::string raw_line;
  std::optional<ConditionRecord> current;
  int header_state = 0;  // 1 = expecting table+path, 2 = expecting timestamps, 3 = in body

  auto flush_current = [&]() {
    if (current.has_value() && !current->table.empty() && !current->field.empty()) {
      records_.push_back(std::move(*current));
      ++added;
    }
    current.reset();
  };

  while (std::getline(in, raw_line)) {
    const std::string line = trim(raw_line);
    if (line.empty() || line[0] == '!') {
      continue;
    }

    if (line.rfind("*DBRECORD", 0) == 0) {
      flush_current();
      current.emplace();
      header_state = 1;
      continue;
    }

    if (!current.has_value()) {
      continue;
    }

    if (header_state == 1) {
      std::istringstream iss(line);
      std::string table;
      std::string full_path;
      iss >> table >> full_path;
      current->table = table;
      if (full_path.size() >= 2 &&
          full_path.compare(full_path.size() - 2, 2, ".W") == 0) {
        current->is_word_update = true;
        current->path = full_path.substr(0, full_path.size() - 2);
      } else {
        current->is_word_update = false;
        current->path = full_path;
      }
      header_state = 2;
      continue;
    }

    if (header_state == 2) {
      const auto parts = split_csv(line);
      if (parts.size() >= 2) {
        current->valid_from.yymmdd = std::stoi(parts[0]);
        current->valid_from.hhmmss = std::stoi(parts[1]);
      }
      if (parts.size() >= 4) {
        current->entered_on.yymmdd = std::stoi(parts[2]);
        current->entered_on.hhmmss = std::stoi(parts[3]);
      }
      header_state = 3;
      continue;
    }

    if (line.rfind("*DBFIELD", 0) == 0) {
      std::istringstream iss(line);
      std::string tag;
      std::string field_name;
      iss >> tag >> field_name;
      current->field = field_name;
      continue;
    }

    // Word update payload line: e.g. "3,SLOW,14., 1.10208"
    const auto parts = split_csv(line);
    if (current->is_word_update && parts.size() >= 4) {
      // parts[2] is 1-based word index (often written with trailing '.', e.g. "14.")
      const std::int32_t word_idx = static_cast<std::int32_t>(std::stof(parts[2]));
      const float val = std::stof(parts[3]);
      current->words[word_idx] = val;
    }
  }

  flush_current();
  return added;
}

std::size_t ConditionDatabase::load_kapack(std::istream& in,
                                           std::size_t max_entries) {
  KapackReader reader(in);
  const auto leaves = reader.list_leaf_entries(max_entries);
  std::size_t added = 0;

  auto is_plausible_yymmdd = [](std::int32_t d) noexcept {
    return d >= 101 && d <= 991231;
  };
  auto is_plausible_hhmmss = [](std::int32_t t) noexcept {
    return t >= 0 && t <= 240000;
  };

  for (const auto& entry : leaves) {
    const std::string key =
        trim(unpack_kapack_key(entry.key_words, reader.header().max_pack));
    if (key.empty() || entry.payload_words.empty()) {
      continue;
    }

    std::vector<std::string> parts;
    std::string cur;
    for (const char ch : key) {
      if (ch == ':' || ch == '/' || std::isspace(static_cast<unsigned char>(ch))) {
        if (!cur.empty()) {
          parts.push_back(std::move(cur));
          cur.clear();
        }
      } else {
        cur.push_back(ch);
      }
    }
    if (!cur.empty()) {
      parts.push_back(std::move(cur));
    }
    if (parts.empty()) {
      continue;
    }

    ConditionRecord rec;
    if (parts.size() >= 3) {
      rec.table = parts[0];
      rec.path = parts[1];
      rec.field = parts[2];
    } else if (parts.size() == 2) {
      rec.table = parts[0];
      rec.path = parts[1];
      rec.field = "DATA";
    } else {
      rec.table = "CALB";
      rec.path = parts[0];
      rec.field = "DATA";
    }

    if (rec.path.size() >= 2 &&
        rec.path.compare(rec.path.size() - 2, 2, ".W") == 0) {
      rec.is_word_update = true;
      rec.path = rec.path.substr(0, rec.path.size() - 2);
    }

    std::size_t data_offset = 0;
    const auto& pw = entry.payload_words;
    if (pw.size() > 4 && is_plausible_yymmdd(pw[0]) &&
        is_plausible_hhmmss(pw[1]) && is_plausible_yymmdd(pw[2]) &&
        is_plausible_hhmmss(pw[3])) {
      rec.valid_from = DbTimestamp{pw[0], pw[1]};
      rec.entered_on = DbTimestamp{pw[2], pw[3]};
      data_offset = 4;
    } else {
      rec.valid_from = DbTimestamp{800101, 0};
      rec.entered_on = DbTimestamp{940101, 0};
    }

    std::int32_t word_idx = 1;
    for (std::size_t i = data_offset; i < pw.size(); ++i, ++word_idx) {
      const float val = std::bit_cast<float>(pw[i]);
      if (std::isfinite(val)) {
        rec.words[word_idx] = val;
      }
    }
    if (!rec.words.empty()) {
      records_.push_back(std::move(rec));
      ++added;
    }
  }
  return added;
}

std::map<std::int32_t, float> ConditionDatabase::resolve_field(
    std::string_view table,
    std::string_view path,
    std::string_view field,
    DbTimestamp query_time,
    DbTimestamp entry_cutoff) const {
  std::vector<const ConditionRecord*> matches;
  for (const auto& rec : records_) {
    if (rec.table != table || rec.path != path || rec.field != field) {
      continue;
    }
    if (rec.valid_from > query_time) {
      continue;
    }
    if (rec.entered_on >= entry_cutoff) {
      continue;
    }
    matches.push_back(&rec);
  }

  std::stable_sort(matches.begin(), matches.end(),
                   [](const ConditionRecord* a, const ConditionRecord* b) {
                     if (a->valid_from != b->valid_from) {
                       return a->valid_from < b->valid_from;
                     }
                     return a->entered_on < b->entered_on;
                   });

  std::map<std::int32_t, float> merged;
  for (const ConditionRecord* rec : matches) {
    if (!rec->is_word_update) {
      merged = rec->words;
    } else {
      for (const auto& [idx, val] : rec->words) {
        merged[idx] = val;
      }
    }
  }
  return merged;
}

std::optional<float> ConditionDatabase::resolve_word(
    std::string_view table,
    std::string_view path,
    std::string_view field,
    std::int32_t one_based_word,
    DbTimestamp query_time,
    DbTimestamp entry_cutoff) const {
  const auto merged = resolve_field(table, path, field, query_time, entry_cutoff);
  const auto it = merged.find(one_based_word);
  if (it == merged.end()) {
    return std::nullopt;
  }
  return it->second;
}

}  // namespace delphipp::db
