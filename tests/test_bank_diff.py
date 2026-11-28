#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit tests for tools/bank_diff.py."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from bank_diff import TolerancePolicy, compare_banks, extract_events  # noqa: E402


class BankDiffTests(unittest.TestCase):
    def _sample_event(self, pz: float = 45.6, dfix_date: int = 20261124) -> list[dict]:
        return [
            {
                "run": 47014,
                "event": 1,
                "banks": [
                    {"idh": "DFIX", "idn": 1, "iwords": {"1": dfix_date}, "fwords": {}},
                    {"idh": "TK", "idn": 7, "iwords": {"1": 13}, "fwords": {"2": pz}},
                ],
            }
        ]

    def test_exact_match_passes(self) -> None:
        errs = compare_banks(self._sample_event(), self._sample_event(), TolerancePolicy())
        self.assertEqual(errs, [])

    def test_float_within_tolerance_passes(self) -> None:
        errs = compare_banks(
            self._sample_event(pz=45.60000),
            self._sample_event(pz=45.60001),
            TolerancePolicy(default_atol=1e-4, default_rtol=1e-4),
        )
        self.assertEqual(errs, [])

    def test_float_outside_tolerance_fails(self) -> None:
        errs = compare_banks(
            self._sample_event(pz=45.6),
            self._sample_event(pz=46.1),
            TolerancePolicy(default_atol=1e-5, default_rtol=1e-4),
        )
        self.assertEqual(len(errs), 1)
        self.assertIn("Bank TK float word 2", errs[0])

    def test_dfix_date_mask_ignored(self) -> None:
        policy = TolerancePolicy(bank_overrides={"DFIX": {"ignore_words": [1]}})
        errs = compare_banks(
            self._sample_event(dfix_date=19990101),
            self._sample_event(dfix_date=20261124),
            policy,
        )
        self.assertEqual(errs, [])

    def test_extract_events_handles_list_and_dict(self) -> None:
        evs = self._sample_event()
        self.assertEqual(extract_events({"events": evs}), evs)
        multi = [{"events": evs}, {"events": evs}]
        self.assertEqual(len(extract_events(multi)), 2)


if __name__ == "__main__":
    unittest.main()
