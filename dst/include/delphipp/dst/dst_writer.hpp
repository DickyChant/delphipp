// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Analysis-ready columnar & HEP event writers for DELPHI++ Short-DST events:
//   - CSV particle table (directly readable by pandas / awkward / pyarrow)
//   - NDJSON (newline-delimited JSON, one DstEvent per line)
//   - HepMC3 Asciiv3 (directly readable by Rivet / HepMC3 tools)

#pragma once

#include <ostream>
#include <span>

#include "delphipp/dst/short_dst.hpp"

namespace delphipp::dst {

void write_csv_header(std::ostream& out);
void write_csv_event(std::ostream& out, const DstEvent& ev);

void write_ndjson_event(std::ostream& out, const DstEvent& ev);

void write_hepmc3_header(std::ostream& out);
void write_hepmc3_event(std::ostream& out, const DstEvent& ev);
void write_hepmc3_footer(std::ostream& out);

}  // namespace delphipp::dst
