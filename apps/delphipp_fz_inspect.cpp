// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// CLI utility to inspect real CERNLIB ZEBRA FZ Exchange files (.fz, .fadsim,
// .fadana, .sdst) produced by DELSIM / delphi-sim-pipeline.

#include <fstream>
#include <iostream>
#include <string>

#include "delphipp/io/fz_scanner.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: delphipp_fz_inspect <file.fz|simana.fadsim|simana.sdst>\n";
    return 2;
  }
  const std::string path = argv[1];
  std::ifstream ifs(path, std::ios::binary);
  if (!ifs) {
    std::cerr << "ERROR: cannot open " << path << '\n';
    return 1;
  }

  const auto s = delphipp::io::scan_fz_stream(ifs);
  std::cout << "{\n"
            << "  \"file\": \"" << path << "\",\n"
            << "  \"file_bytes\": " << s.file_bytes << ",\n"
            << "  \"detected_zebra_fz\": " << (s.detected_zebra_fz ? "true" : "false") << ",\n"
            << "  \"endianness\": \"" << (s.is_big_endian ? "big" : "little") << "\",\n"
            << "  \"physical_record_words\": " << s.physical_record_words << ",\n"
            << "  \"steering_blocks\": " << s.steering_blocks << ",\n"
            << "  \"fast_blocks\": " << s.fast_blocks << ",\n"
            << "  \"logical_records\": " << s.logical_records << ",\n"
            << "  \"bank_id_hits\": {";
  bool first = true;
  for (const auto& [k, v] : s.bank_id_counts) {
    if (!first) std::cout << ", ";
    first = false;
    std::cout << "\"" << k << "\": " << v;
  }
  std::cout << "}\n}\n";
  return s.detected_zebra_fz ? 0 : 1;
}
