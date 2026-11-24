// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/io/fz_exchange.hpp"

#include <cstring>
#include <stdexcept>

namespace delphipp::io {

std::uint32_t read_be32(std::span<const std::byte, 4> b) noexcept {
  return (static_cast<std::uint32_t>(std::to_integer<std::uint8_t>(b[0])) << 24) |
         (static_cast<std::uint32_t>(std::to_integer<std::uint8_t>(b[1])) << 16) |
         (static_cast<std::uint32_t>(std::to_integer<std::uint8_t>(b[2])) << 8) |
         (static_cast<std::uint32_t>(std::to_integer<std::uint8_t>(b[3])));
}

void write_be32(std::uint32_t v, std::span<std::byte, 4> out) noexcept {
  out[0] = static_cast<std::byte>((v >> 24) & 0xFFu);
  out[1] = static_cast<std::byte>((v >> 16) & 0xFFu);
  out[2] = static_cast<std::byte>((v >> 8) & 0xFFu);
  out[3] = static_cast<std::byte>(v & 0xFFu);
}

FzExchangeWriter::FzExchangeWriter(std::ostream& out,
                                   std::size_t physical_record_words)
    : out_(out), phys_words_(physical_record_words) {
  if (phys_words_ <= kSteeringHeaderWords + 4) {
    throw std::invalid_argument("FZ physical record size too small");
  }
  start_new_block();
}

void FzExchangeWriter::start_new_block() {
  block_.clear();
  block_.reserve(phys_words_);
  for (const auto w : kFzSyncWords) {
    block_.push_back(w);
  }
  block_.push_back(static_cast<std::uint32_t>(phys_words_));
  block_.push_back(phys_counter_);
  block_.push_back(0u);  // steering flags
  block_.push_back(static_cast<std::uint32_t>(kSteeringHeaderWords));  // pilot offset
}

void FzExchangeWriter::flush_block() {
  if (block_.size() <= kSteeringHeaderWords) {
    return;
  }
  while (block_.size() < phys_words_) {
    block_.push_back(0u);
  }
  std::vector<std::byte> raw(phys_words_ * 4);
  for (std::size_t i = 0; i < phys_words_; ++i) {
    write_be32(block_[i], std::span<std::byte, 4>(raw.data() + i * 4, 4));
  }
  out_.write(reinterpret_cast<const char*>(raw.data()),
             static_cast<std::streamsize>(raw.size()));
  ++phys_counter_;
  start_new_block();
}

void FzExchangeWriter::append_word(std::uint32_t word) {
  if (block_.size() >= phys_words_) {
    flush_block();
  }
  block_.push_back(word);
}

void FzExchangeWriter::write_record(const FzLogicalRecord& rec) {
  const std::uint32_t nuh = static_cast<std::uint32_t>(rec.user_header.size());
  const std::uint32_t npay = static_cast<std::uint32_t>(rec.payload.size());
  // Pilot header: [total_words_following, record_type, nuh, npay]
  const std::uint32_t total_len = 2u + nuh + npay;
  append_word(total_len);
  append_word(static_cast<std::uint32_t>(rec.type));
  append_word(nuh);
  append_word(npay);
  for (const std::int32_t w : rec.user_header) {
    append_word(static_cast<std::uint32_t>(w));
  }
  for (const std::int32_t w : rec.payload) {
    append_word(static_cast<std::uint32_t>(w));
  }
}

void FzExchangeWriter::finish() {
  write_record(FzLogicalRecord{.type = FzRecordType::EndOfFile});
  flush_block();
}

FzExchangeReader::FzExchangeReader(std::istream& in) : in_(in) {}

bool FzExchangeReader::load_next_block() {
  std::array<std::byte, kSteeringHeaderWords * 4> hdr_bytes{};
  if (!in_.read(reinterpret_cast<char*>(hdr_bytes.data()),
                static_cast<std::streamsize>(hdr_bytes.size()))) {
    return false;
  }

  std::array<std::uint32_t, kSteeringHeaderWords> hdr{};
  for (std::size_t i = 0; i < kSteeringHeaderWords; ++i) {
    hdr[i] = read_be32(std::span<const std::byte, 4>(hdr_bytes.data() + i * 4, 4));
  }

  for (std::size_t i = 0; i < kFzSyncWords.size(); ++i) {
    if (hdr[i] != kFzSyncWords[i]) {
      throw std::runtime_error("FZ Exchange sync word mismatch");
    }
  }

  const std::size_t phys_words = hdr[4];
  if (phys_words <= kSteeringHeaderWords) {
    throw std::runtime_error("Invalid FZ physical record length");
  }

  const std::size_t rem_words = phys_words - kSteeringHeaderWords;
  std::vector<std::byte> body_bytes(rem_words * 4);
  if (!in_.read(reinterpret_cast<char*>(body_bytes.data()),
                static_cast<std::streamsize>(body_bytes.size()))) {
    throw std::runtime_error("Truncated FZ physical record");
  }

  block_.resize(phys_words);
  for (std::size_t i = 0; i < kSteeringHeaderWords; ++i) {
    block_[i] = hdr[i];
  }
  for (std::size_t i = 0; i < rem_words; ++i) {
    block_[kSteeringHeaderWords + i] =
        read_be32(std::span<const std::byte, 4>(body_bytes.data() + i * 4, 4));
  }

  cursor_ = kSteeringHeaderWords;
  ++phys_read_;
  return true;
}

std::optional<std::uint32_t> FzExchangeReader::read_word() {
  if (cursor_ >= block_.size()) {
    if (!load_next_block()) {
      return std::nullopt;
    }
  }
  return block_[cursor_++];
}

std::optional<FzLogicalRecord> FzExchangeReader::next_record() {
  while (true) {
    const auto len_opt = read_word();
    if (!len_opt.has_value()) {
      return std::nullopt;
    }
    const std::uint32_t total_len = *len_opt;
    if (total_len == 0u) {
      // Zero padding to end of current physical block
      cursor_ = block_.size();
      continue;
    }

    const auto type_opt = read_word();
    const auto nuh_opt = read_word();
    const auto npay_opt = read_word();
    if (!type_opt || !nuh_opt || !npay_opt) {
      throw std::runtime_error("Truncated FZ logical record pilot");
    }

    const auto rtype = static_cast<FzRecordType>(*type_opt);
    const std::uint32_t nuh = *nuh_opt;
    const std::uint32_t npay = *npay_opt;

    if (total_len != 2u + nuh + npay) {
      throw std::runtime_error("FZ logical record length mismatch");
    }

    FzLogicalRecord rec;
    rec.type = rtype;
    rec.user_header.reserve(nuh);
    for (std::uint32_t i = 0; i < nuh; ++i) {
      const auto w = read_word();
      if (!w) {
        throw std::runtime_error("Truncated FZ user header");
      }
      rec.user_header.push_back(static_cast<std::int32_t>(*w));
    }

    rec.payload.reserve(npay);
    for (std::uint32_t i = 0; i < npay; ++i) {
      const auto w = read_word();
      if (!w) {
        throw std::runtime_error("Truncated FZ payload");
      }
      rec.payload.push_back(static_cast<std::int32_t>(*w));
    }

    if (rtype == FzRecordType::Padding) {
      continue;
    }
    if (rtype == FzRecordType::EndOfFile) {
      return std::nullopt;
    }
    return rec;
  }
}

}  // namespace delphipp::io
