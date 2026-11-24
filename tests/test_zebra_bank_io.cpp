// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>

#include "delphipp/io/fz_exchange.hpp"
#include "delphipp/io/zebra_bank_io.hpp"
#include "delphipp/zebra_store.hpp"

int main() {
  using namespace delphipp;
  using namespace delphipp::io;

  OwnedZebraStore src_store(512);
  auto src = src_store.view();

  const std::int32_t l_src = 100;
  src.set_iq(l_src - 5, 42);
  std::int32_t idh = 0;
  std::memcpy(&idh, "EGRA", 4);
  src.set_iq(l_src - 4, idh);
  src.set_iq(l_src - 3, 2);   // nl
  src.set_iq(l_src - 2, 1);   // ns
  src.set_iq(l_src - 1, 3);   // nd
  src.set_iq(l_src, 0x05);    // status
  src.set_lq(l_src - 1, 200); // link 1
  src.set_lq(l_src - 2, 300); // link 2
  src.set_iq(l_src + 1, 29);
  src.set_q(l_src + 2, 91.1876f);
  src.set_iq(l_src + 3, -7);

  std::stringstream fz(std::ios::in | std::ios::out | std::ios::binary);
  {
    FzExchangeWriter writer(fz, 64);
    writer.write_record(pack_zebra_bank(src, l_src, {47014, 1}));
    writer.finish();
  }

  fz.seekg(0);
  FzExchangeReader reader(fz);
  const auto rec = reader.next_record();
  assert(rec.has_value());

  OwnedZebraStore dst_store(512);
  auto dst = dst_store.view();
  const std::int32_t l_dst = unpack_zebra_bank(*rec, dst, 150);

  assert(dst.bank_idn(l_dst) == 42);
  assert(dst.bank_idh(l_dst) == "EGRA");
  assert(dst.bank_nl(l_dst) == 2);
  assert(dst.bank_ns(l_dst) == 1);
  assert(dst.bank_nd(l_dst) == 3);
  assert(dst.bank_status(l_dst) == 0x05);
  assert(dst.bank_link(l_dst, 1) == 200);
  assert(dst.bank_link(l_dst, 2) == 300);
  assert(dst.bank_idata(l_dst, 1) == 29);
  assert(std::fabs(dst.bank_fdata(l_dst, 2) - 91.1876f) < 1e-5f);
  assert(dst.bank_idata(l_dst, 3) == -7);

  std::cout << "test_zebra_bank_io: PASS\n";
  return 0;
}
