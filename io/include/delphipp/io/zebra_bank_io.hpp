// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Serialize / materialize ZEBRA bank structures between ZebraStoreView and
// FzLogicalRecord payloads.

#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "delphipp/io/fz_exchange.hpp"
#include "delphipp/zebra_store.hpp"

namespace delphipp::io {

// Pack a single ZEBRA bank at 1-based address `l` from `store` into an
// FzLogicalRecord payload (preserving nl, ns, nd, idn, idh, status, links, data).
[[nodiscard]] FzLogicalRecord pack_zebra_bank(const ZebraStoreView& store,
                                              std::int32_t l,
                                              std::vector<std::int32_t> user_header = {});

// Materialize a packed ZEBRA bank from `rec` into `store` starting at 1-based
// address `target_l` (requiring target_l > nl + 5 and sufficient store capacity).
// Returns `target_l`.
std::int32_t unpack_zebra_bank(const FzLogicalRecord& rec,
                               ZebraStoreView& store,
                               std::int32_t target_l);

}  // namespace delphipp::io
