// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Pure C++20 CERNLIB ZEBRA FZ Exchange binary format reader and writer.
// Implements 32-bit big-endian physical steering blocks (sync words
// 0x0123CDEF, 0x80708070, 0x4321ABCD, 0x80618061), multi-block logical
// records, user headers, and payload word extraction.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <optional>
#include <ostream>
#include <span>
#include <vector>

namespace delphipp::io {

// Standard CERNLIB ZEBRA FZ Exchange steering block sync words
inline constexpr std::array<std::uint32_t, 4> kFzSyncWords = {
    0x0123CDEFu,
    0x80708070u,
    0x4321ABCDu,
    0x80618061u,
};

inline constexpr std::size_t kDefaultPhysicalRecordWords = 900;  // 3600 bytes
inline constexpr std::size_t kSteeringHeaderWords = 8;

enum class FzRecordType : std::uint32_t {
  StartOfRun = 1,
  EventNormal = 2,
  EventContinue = 3,
  Padding = 4,
  EndOfFile = 5,
  EndOfRun = 6,
};

struct FzLogicalRecord {
  FzRecordType type = FzRecordType::EventNormal;
  std::vector<std::int32_t> user_header;
  std::vector<std::int32_t> payload;
};

// Big-endian 32-bit word serialization helpers
[[nodiscard]] std::uint32_t read_be32(std::span<const std::byte, 4> bytes) noexcept;
void write_be32(std::uint32_t value, std::span<std::byte, 4> out) noexcept;

class FzExchangeWriter {
 public:
  explicit FzExchangeWriter(
      std::ostream& out,
      std::size_t physical_record_words = kDefaultPhysicalRecordWords);

  void write_record(const FzLogicalRecord& rec);
  void finish();  // Flushes pending physical block and writes EndOfFile marker

 private:
  std::ostream& out_;
  std::size_t phys_words_;
  std::uint32_t phys_counter_ = 0;
  std::vector<std::uint32_t> block_;

  void start_new_block();
  void flush_block();
  void append_word(std::uint32_t word);
};

class FzExchangeReader {
 public:
  explicit FzExchangeReader(std::istream& in);

  // Reads the next non-padding logical record, or std::nullopt at EOF / EndOfFile.
  [[nodiscard]] std::optional<FzLogicalRecord> next_record();

  [[nodiscard]] std::uint32_t physical_records_read() const noexcept {
    return phys_read_;
  }

 private:
  std::istream& in_;
  std::vector<std::uint32_t> block_;
  std::size_t cursor_ = 0;
  std::uint32_t phys_read_ = 0;

  bool load_next_block();
  std::optional<std::uint32_t> read_word();
};

}  // namespace delphipp::io
