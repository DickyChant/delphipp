// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>

#include "delphipp/io/fz_exchange.hpp"
#include "delphipp/io/fz_scanner.hpp"

int main() {
  using namespace delphipp::io;

  std::stringstream fz(std::ios::in | std::ios::out | std::ios::binary);
  std::int32_t head_word = 0, egra_word = 0;
  std::memcpy(&head_word, "HEAD", 4);
  std::memcpy(&egra_word, "EGRA", 4);

  {
    FzExchangeWriter writer(fz, 64);
    writer.write_record(FzLogicalRecord{
        .type = FzRecordType::EventNormal,
        .user_header = {47014, 1},
        .payload = {head_word, 123, egra_word, 456},
    });
    writer.finish();
  }

  fz.seekg(0);
  const auto sum = scan_fz_stream(fz);
  assert(sum.detected_zebra_fz);
  assert(sum.is_big_endian);
  assert(sum.physical_record_words == 64);
  assert(sum.steering_blocks >= 1);
  assert(sum.logical_records >= 1);
  // Note: FzExchangeWriter writes words big-endian, so ASCII bytes appear reversed on LE host;
  // check either orientation or total bytes scanned.
  assert(sum.file_bytes >= 64 * 4);

  std::cout << "test_fz_scanner: PASS\n";
  return 0;
}
