// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "delphipp/io/zebra_bank_io.hpp"

#include <stdexcept>

namespace delphipp::io {

FzLogicalRecord pack_zebra_bank(const ZebraStoreView& store,
                                std::int32_t l,
                                std::vector<std::int32_t> user_header) {
  const std::int32_t idn = store.bank_idn(l);
  const std::int32_t idh = store.iq(l - 4);
  const std::int32_t nl = store.bank_nl(l);
  const std::int32_t ns = store.bank_ns(l);
  const std::int32_t nd = store.bank_nd(l);
  const std::int32_t status = store.bank_status(l);

  FzLogicalRecord rec;
  rec.type = FzRecordType::EventNormal;
  rec.user_header = std::move(user_header);
  rec.payload.reserve(static_cast<std::size_t>(6 + nl + nd));
  rec.payload.push_back(idn);
  rec.payload.push_back(idh);
  rec.payload.push_back(nl);
  rec.payload.push_back(ns);
  rec.payload.push_back(nd);
  rec.payload.push_back(status);
  for (std::int32_t j = 1; j <= nl; ++j) {
    rec.payload.push_back(store.bank_link(l, j));
  }
  for (std::int32_t k = 1; k <= nd; ++k) {
    rec.payload.push_back(store.bank_idata(l, k));
  }
  return rec;
}

std::int32_t unpack_zebra_bank(const FzLogicalRecord& rec,
                               ZebraStoreView& store,
                               std::int32_t target_l) {
  if (rec.payload.size() < 6) {
    throw std::runtime_error("ZEBRA bank payload shorter than 6 header words");
  }
  const std::int32_t idn = rec.payload[0];
  const std::int32_t idh = rec.payload[1];
  const std::int32_t nl = rec.payload[2];
  const std::int32_t ns = rec.payload[3];
  const std::int32_t nd = rec.payload[4];
  const std::int32_t status = rec.payload[5];

  if (nl < 0 || nd < 0 ||
      rec.payload.size() != static_cast<std::size_t>(6 + nl + nd)) {
    throw std::runtime_error("ZEBRA bank payload length mismatch");
  }

  store.set_lq(target_l, 0);
  store.set_lq(target_l + 1, 0);
  store.set_iq(target_l - 5, idn);
  store.set_iq(target_l - 4, idh);
  store.set_iq(target_l - 3, nl);
  store.set_iq(target_l - 2, ns);
  store.set_iq(target_l - 1, nd);
  store.set_iq(target_l, status);

  std::size_t idx = 6;
  for (std::int32_t j = 1; j <= nl; ++j) {
    store.set_lq(target_l - j, rec.payload[idx++]);
  }
  for (std::int32_t k = 1; k <= nd; ++k) {
    store.set_iq(target_l + k, rec.payload[idx++]);
  }
  return target_l;
}

}  // namespace delphipp::io
