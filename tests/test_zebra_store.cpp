// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>

#include "delphipp/zebra_store.hpp"

int main() {
  delphipp::OwnedZebraStore store(1024);
  auto view = store.view();

  // Verify EQUIVALENCE (Q(1), IQ(1), LQ(9))
  view.set_iq(1, 0x12345678);
  assert(view.lq(9) == 0x12345678);

  // Verify Q(k) bit-cast equivalence with IQ(k)
  view.set_q(42, 91.1876f);
  assert(std::fabs(view.q(42) - 91.1876f) < 1e-6f);

  // Construct a minimal synthetic ZEBRA bank at L = 100
  const std::int32_t l = 100;
  view.set_lq(l, 0);      // next
  view.set_lq(l + 1, 50); // up
  view.set_iq(l - 5, 7);  // idn
  std::int32_t idh = 0;
  std::memcpy(&idh, "TK  ", 4);
  view.set_iq(l - 4, idh); // idh
  view.set_iq(l - 3, 2);   // nl
  view.set_iq(l - 2, 1);   // ns
  view.set_iq(l - 1, 3);   // nd
  view.set_lq(l - 1, 200); // structural link 1
  view.set_iq(l + 1, 13);  // data word 1 (int)
  view.set_q(l + 2, 45.6f); // data word 2 (float)

  assert(view.bank_up(l) == 50);
  assert(view.bank_idn(l) == 7);
  assert(view.bank_idh(l) == "TK  ");
  assert(view.bank_nl(l) == 2);
  assert(view.bank_ns(l) == 1);
  assert(view.bank_nd(l) == 3);
  assert(view.bank_link(l, 1) == 200);
  assert(view.bank_idata(l, 1) == 13);
  assert(std::fabs(view.bank_fdata(l, 2) - 45.6f) < 1e-6f);

  std::cout << "test_zebra_store: PASS\n";
  return 0;
}
