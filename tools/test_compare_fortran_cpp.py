#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest

from tools.compare_fortran_cpp import compare_reports, parse_fortran_dump


SAMPLE_SDST = """\
=== FORTRAN CERNLIB FZ DUMP: simana.sdst
  EVT #    1 run=  -47014 evt=       1 top=DTOP nd=     35 banks=    12
  EVT #    2 run=  -47014 evt=       2 top=DTOP nd=     35 banks=    14
Total events read by CERNLIB FZIN:      2 LREC= 3840
"""

SAMPLE_FADSIM_WITH_RUN = """\
=== FORTRAN CERNLIB FZ DUMP: simana.fadsim
  EVT #    1 run=  -47014 evt=       0 top=RUN  nd=     71 banks=     1
  EVT #    2 run=  -47014 evt=       1 top=TOP  nd=      1 banks=    18
  EVT #    3 run=  -47014 evt=       2 top=TOP  nd=      1 banks=    19
Total events read by CERNLIB FZIN:      3 LREC= 3840
"""


class CompareFortranCppTest(unittest.TestCase):
    def test_sdst_matching_reports(self) -> None:
        fevs, ftotal = parse_fortran_dump(SAMPLE_SDST)
        self.assertEqual(ftotal, 2)
        self.assertEqual(len(fevs), 2)
        self.assertEqual(fevs[0].run, -47014)
        self.assertEqual(fevs[0].top, "DTOP")

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
        ok, report = compare_reports(SAMPLE_SDST, cpp_summary, "simana.sdst")
        self.assertTrue(ok, report)

    def test_fadsim_run_preamble_excluded_from_physics_count(self) -> None:
        cpp_summary = {
            "physical_record_words": 3840,
            "events_count": 2,
            "events": [
                {"index": 1, "run": 47014, "event": 1, "zebra_banks": 18},
                {"index": 2, "run": 47014, "event": 2, "zebra_banks": 19},
            ],
        }
        ok, report = compare_reports(SAMPLE_FADSIM_WITH_RUN, cpp_summary, "simana.fadsim")
        self.assertTrue(ok, report)

    def test_event_count_mismatch_fails(self) -> None:
        cpp_summary = {"physical_record_words": 3840, "events_count": 1, "events": []}
        ok, _ = compare_reports(SAMPLE_SDST, cpp_summary, "bad")
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main()
