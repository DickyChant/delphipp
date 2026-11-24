#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest

from tools.compare_fortran_cpp import compare_reports, parse_fortran_dump


SAMPLE_FORTRAN = """\
=== FORTRAN CERNLIB FZ DUMP: simana.sdst
  EVT #   1 run= 47014 evt=     1 top=HEAD nd=   24 banks=    12
  EVT #   2 run= 47014 evt=     2 top=HEAD nd=   24 banks=    14
Total events read by CERNLIB FZIN:      2 LREC= 3840
"""


class CompareFortranCppTest(unittest.TestCase):
    def test_matching_reports(self) -> None:
        fevs, ftotal = parse_fortran_dump(SAMPLE_FORTRAN)
        self.assertEqual(ftotal, 2)
        self.assertEqual(len(fevs), 2)
        self.assertEqual(fevs[0].run, 47014)

        cpp_summary = {
            "physical_record_words": 3840,
            "events_count": 2,
            "events": [
                {
                    "index": 1,
                    "run": 47014,
                    "event": 1,
                    "zebra_banks": 15,
                    "raw_blocklets": 0,
                    "dst_particles": 10,
                    "visible_energy_gev": 85.4,
                },
                {
                    "index": 2,
                    "run": 47014,
                    "event": 2,
                    "zebra_banks": 16,
                    "raw_blocklets": 0,
                    "dst_particles": 11,
                    "visible_energy_gev": 90.1,
                },
            ],
        }
        ok, report = compare_reports(SAMPLE_FORTRAN, cpp_summary, "simana.sdst")
        self.assertTrue(ok, report)

    def test_event_count_mismatch_fails(self) -> None:
        cpp_summary = {"physical_record_words": 3840, "events_count": 1, "events": []}
        ok, _ = compare_reports(SAMPLE_FORTRAN, cpp_summary, "bad")
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main()
