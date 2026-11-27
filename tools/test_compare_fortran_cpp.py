#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import unittest

from tools.compare_fortran_cpp import compare_reports, parse_fortran_dump


SAMPLE_SDST = """\
=== FORTRAN CERNLIB FZ DUMP: simana.sdst
  EVT #    1 run=  -47014 evt=       1 top=DTOP nd=     35 banks=    12
    BANK_STAT npa=   10 nchg=    6 nneu=    4 sum_e=    85.400 sum_p=    85.100
  EVT #    2 run=  -47014 evt=       2 top=DTOP nd=     35 banks=    14
    BANK_STAT npa=   11 nchg=    7 nneu=    4 sum_e=    90.100 sum_p=    89.700
Total events read by CERNLIB FZIN:      2 LREC= 3840
"""

SAMPLE_FADSIM_MIXED = """\
=== FORTRAN CERNLIB FZ DUMP: simana.fadsim
  EVT #    1 run=  -47014 evt=       0 top=RUN  nd=     71 banks=     1
  EVT #    2 run=  -47014 evt=       0 top=SIM_ nd=      1 banks=     4
  EVT #    3 run=  -47014 evt=       0 top=TOP  nd=      1 banks=     2
    RAW_STAT nrawd=     0 npmb=   0
  EVT #    4 run=  -47014 evt=       1 top=SIM_ nd=      1 banks=    18
  EVT #    5 run=  -47014 evt=       1 top=TOP  nd=      1 banks=    19
    RAW_STAT nrawd=    14 npmb=  12
Total events read by CERNLIB FZIN:      5 LREC= 3840
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
                    "run": -47014,
                    "event": 1,
                    "top_bank": "DTOP",
                    "zebra_banks": 15,
                    "raw_blocklets": 0,
                    "dst_particles": 10,
                    "visible_energy_gev": 85.4,
                },
                {
                    "index": 2,
                    "run": -47014,
                    "event": 2,
                    "top_bank": "DTOP",
                    "zebra_banks": 16,
                    "raw_blocklets": 0,
                    "dst_particles": 11,
                    "visible_energy_gev": 90.1,
                },
            ],
        }
        ok, report = compare_reports(SAMPLE_SDST, cpp_summary, "simana.sdst")
        self.assertTrue(ok, report)

    def test_fadsim_mixed_sim_and_top_with_evt0(self) -> None:
        cpp_summary = {
            "physical_record_words": 3840,
            "events_count": 4,
            "events": [
                {"index": 1, "run": -47014, "event": 0, "top_bank": "SIM_", "zebra_banks": 69, "raw_blocklets": 0},
                {"index": 2, "run": -47014, "event": 0, "top_bank": "TOP", "zebra_banks": 4, "raw_blocklets": 0},
                {"index": 3, "run": -47014, "event": 1, "top_bank": "SIM_", "zebra_banks": 109, "raw_blocklets": 0},
                {"index": 4, "run": -47014, "event": 1, "top_bank": "TOP", "zebra_banks": 1132, "raw_blocklets": 14},
            ],
        }
        ok, report = compare_reports(SAMPLE_FADSIM_MIXED, cpp_summary, "simana.fadsim")
        self.assertTrue(ok, report)

    def test_bank_and_raw_stat_mismatch_fails(self) -> None:
        bad_sdst = {
            "physical_record_words": 3840,
            "events_count": 2,
            "events": [
                {"index": 1, "run": -47014, "event": 1, "top_bank": "DTOP", "dst_particles": 9, "visible_energy_gev": 85.4},
                {"index": 2, "run": -47014, "event": 2, "top_bank": "DTOP", "dst_particles": 11, "visible_energy_gev": 90.1},
            ],
        }
        ok_sdst, _ = compare_reports(SAMPLE_SDST, bad_sdst, "bad_npa")
        self.assertFalse(ok_sdst)

        bad_raw = {
            "physical_record_words": 3840,
            "events_count": 4,
            "events": [
                {"index": 1, "run": -47014, "event": 0, "top_bank": "SIM_", "raw_blocklets": 0},
                {"index": 2, "run": -47014, "event": 0, "top_bank": "TOP", "raw_blocklets": 0},
                {"index": 3, "run": -47014, "event": 1, "top_bank": "SIM_", "raw_blocklets": 0},
                {"index": 4, "run": -47014, "event": 1, "top_bank": "TOP", "raw_blocklets": 99},
            ],
        }
        ok_raw, _ = compare_reports(SAMPLE_FADSIM_MIXED, bad_raw, "bad_nrawd")
        self.assertFalse(ok_raw)

    def test_event_number_mismatch_fails(self) -> None:
        cpp_summary = {
            "physical_record_words": 3840,
            "events_count": 2,
            "events": [
                {"index": 1, "run": -47014, "event": 99, "top_bank": "DTOP"},
                {"index": 2, "run": -47014, "event": 2, "top_bank": "DTOP"},
            ],
        }
        ok, _ = compare_reports(SAMPLE_SDST, cpp_summary, "bad_evt")
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main()
