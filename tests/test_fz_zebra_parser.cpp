// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <bit>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>

#include "delphipp/io/fz_exchange.hpp"
#include "delphipp/io/fz_zebra_parser.hpp"
#include "delphipp/io/zebra_bank_io.hpp"
#include "delphipp/zebra_store.hpp"

namespace {

std::int32_t ascii_to_word(const char s[4]) {
  return (static_cast<std::int32_t>(static_cast<unsigned char>(s[0])) << 24) |
         (static_cast<std::int32_t>(static_cast<unsigned char>(s[1])) << 16) |
         (static_cast<std::int32_t>(static_cast<unsigned char>(s[2])) << 8) |
         (static_cast<std::int32_t>(static_cast<unsigned char>(s[3])));
}

}  // namespace

int main() {
  using namespace delphipp;
  using namespace delphipp::io;

  // 1. Test pack_zebra_bank round-trip via FzExchangeWriter -> parse_fz_stream
  {
    OwnedZebraStore store(2048);
    auto v = store.view();
    const std::int32_t l = 100;
    v.set_iq(l - 5, 1);
    v.set_iq(l - 4, ascii_to_word("HEAD"));
    v.set_iq(l - 3, 2);
    v.set_iq(l - 2, 1);
    v.set_iq(l - 1, 5);
    v.set_iq(l, 0);
    v.set_lq(l - 1, 0);
    v.set_lq(l - 2, 0);
    v.set_iq(l + 1, 100);     // version
    v.set_iq(l + 2, 940526);  // date
    v.set_iq(l + 3, 120000);  // time
    v.set_iq(l + 4, 47014);   // NRUN
    v.set_iq(l + 5, 42);      // IEVT

    std::stringstream fz(std::ios::in | std::ios::out | std::ios::binary);
    {
      FzExchangeWriter w(fz, 64);
      w.write_record(pack_zebra_bank(v, l, {47014, 42}));
      w.finish();
    }
    fz.seekg(0);
    const auto parsed = parse_fz_stream(fz);
    assert(parsed.detected_zebra_fz);
    assert(parsed.is_big_endian);
    assert(parsed.events.size() == 1);
    assert(parsed.events[0].run_number == 47014);
    assert(parsed.events[0].event_number == 42);
    assert(parsed.events[0].banks.size() == 1);
    assert(parsed.events[0].banks[0].idh == "HEAD");
  }

  // 2. Test raw CERNLIB ZEBRA FZ physical steering + fast block stream
  //    containing two events delimited by HEAD banks plus EGRA and PA banks.
  {
    const std::size_t phys_words = 32;
    std::vector<std::uint32_t> raw_words;

    auto emit_steering = [&](std::uint32_t blk_idx) {
      for (const auto sw : kFzSyncWords) {
        raw_words.push_back(sw);
      }
      raw_words.push_back(static_cast<std::uint32_t>(phys_words));
      raw_words.push_back(blk_idx);
      raw_words.push_back(0u);
      raw_words.push_back(8u);
    };

    emit_steering(0);
    // Pilot dummy header (4 words) before ZEBRA material
    raw_words.push_back(40u);
    raw_words.push_back(2u);
    raw_words.push_back(0u);
    raw_words.push_back(0u);

    auto append_bank = [&](const char idh[4], std::int32_t nl, std::int32_t ns,
                           const std::vector<std::int32_t>& data) {
      raw_words.push_back(1u);
      raw_words.push_back(static_cast<std::uint32_t>(ascii_to_word(idh)));
      raw_words.push_back(static_cast<std::uint32_t>(nl));
      raw_words.push_back(static_cast<std::uint32_t>(ns));
      raw_words.push_back(static_cast<std::uint32_t>(data.size()));
      raw_words.push_back(0u);
      for (std::int32_t d : data) {
        raw_words.push_back(static_cast<std::uint32_t>(d));
      }
    };

    // Event 1: HEAD + EGRA
    append_bank("HEAD", 0, 0, {1, 940526, 120000, 47014, 1});
    append_bank("EGRA", 0, 0, {0, (22 << 8) | (1 << 16), 0, 0, 1, 0, 0, 0});
    // Event 2: HEAD + PA
    append_bank("HEAD", 0, 0, {1, 940526, 120001, 47014, 2});
    // Minimal PA bank: 1 module (N+1=1 -> standard block starts at word 6)
    std::vector<std::int32_t> pa_data(32, 0);
    pa_data[0] = 101;
    pa_data[3] = 1;   // N+1 = 1
    pa_data[4] = 26;  // standard length
    pa_data[5 + 2] = std::bit_cast<std::int32_t>(12.5f);   // px
    pa_data[5 + 3] = std::bit_cast<std::int32_t>(-3.2f);   // py
    pa_data[5 + 4] = std::bit_cast<std::int32_t>(18.0f);   // pz
    pa_data[5 + 5] = std::bit_cast<std::int32_t>(22.15f);  // E
    pa_data[5 + 7] = std::bit_cast<std::int32_t>(1.0f);    // charge=+
    append_bank("PA  ", 0, 0, pa_data);

    // Pad out block 1 as fast block if > 32 words
    while (raw_words.size() < phys_words * 2) {
      raw_words.push_back(0u);
    }

    std::string bytes(raw_words.size() * 4, '\0');
    for (std::size_t i = 0; i < raw_words.size(); ++i) {
      const std::uint32_t v = raw_words[i];
      bytes[i * 4 + 0] = static_cast<char>((v >> 24) & 0xFFu);
      bytes[i * 4 + 1] = static_cast<char>((v >> 16) & 0xFFu);
      bytes[i * 4 + 2] = static_cast<char>((v >> 8) & 0xFFu);
      bytes[i * 4 + 3] = static_cast<char>(v & 0xFFu);
    }

    std::istringstream fz(bytes, std::ios::binary);
    const auto parsed = parse_fz_stream(fz);
    assert(parsed.detected_zebra_fz);
    assert(parsed.events.size() == 2);
    assert(parsed.events[0].run_number == 47014);
    assert(parsed.events[0].event_number == 1);
    assert(parsed.events[0].find_first("EGRA") != nullptr);
    assert(parsed.events[1].event_number == 2);
    assert(parsed.events[1].find_first("PA  ") != nullptr);

    OwnedZebraStore zs(4096);
    auto zs_view = zs.view();
    const auto heads = parsed.events[0].populate_zebra_store(zs_view);
    assert(heads.contains("HEAD"));
    assert(heads.contains("EGRA"));
    assert(zs_view.bank_nd(heads.at("EGRA")) == 8);
  }

  std::cout << "test_fz_zebra_parser: PASS\n";
  return 0;
}
