// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <vector>

#include "delphipp/io/fz_exchange.hpp"

int main() {
  using namespace delphipp::io;

  // Use a small physical block (32 words) to exercise multi-block spanning records
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);

  FzLogicalRecord sor{
      .type = FzRecordType::StartOfRun,
      .user_header = {47014, 940526},
      .payload = {100, 200, 300},
  };

  FzLogicalRecord evt{
      .type = FzRecordType::EventNormal,
      .user_header = {47014, 1},
      .payload = {},
  };
  for (std::int32_t i = 0; i < 80; ++i) {
    evt.payload.push_back(i * 17 - 500);
  }

  {
    FzExchangeWriter writer(stream, 32);
    writer.write_record(sor);
    writer.write_record(evt);
    writer.finish();
  }

  stream.seekg(0);
  FzExchangeReader reader(stream);

  const auto r1 = reader.next_record();
  assert(r1.has_value());
  assert(r1->type == FzRecordType::StartOfRun);
  assert(r1->user_header == sor.user_header);
  assert(r1->payload == sor.payload);

  const auto r2 = reader.next_record();
  assert(r2.has_value());
  assert(r2->type == FzRecordType::EventNormal);
  assert(r2->user_header == evt.user_header);
  assert(r2->payload == evt.payload);

  // Next read hits EndOfFile and returns nullopt
  const auto r3 = reader.next_record();
  assert(!r3.has_value());
  // Spanning 80-word payload over 32-word physical blocks requires >= 4 blocks
  assert(reader.physical_records_read() >= 4);

  std::cout << "test_fz_exchange: PASS\n";
  return 0;
}
