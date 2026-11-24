#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
"""Cross-check CERNLIB Fortran FZ reference dump vs DELPHI++ C++20 reco summary."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

EVT_RE = re.compile(
    r"EVT\s+#\s*(\d+)\s+run=\s*(-?\d+)\s+evt=\s*(-?\d+)\s+top=(\S+)\s+nd=\s*(\d+)\s+banks=\s*(\d+)"
)
TOTAL_RE = re.compile(r"Total events read by CERNLIB FZIN:\s*(\d+)")


@dataclass
class FortranEvent:
    index: int
    run: int
    event: int
    top: str
    nd: int
    banks: int


def parse_fortran_dump(text: str) -> tuple[list[FortranEvent], int]:
    events: list[FortranEvent] = []
    for line in text.splitlines():
        m = EVT_RE.search(line)
        if m:
            events.append(
                FortranEvent(
                    index=int(m.group(1)),
                    run=int(m.group(2)),
                    event=int(m.group(3)),
                    top=m.group(4),
                    nd=int(m.group(5)),
                    banks=int(m.group(6)),
                )
            )
    mt = TOTAL_RE.search(text)
    total = int(mt.group(1)) if mt else len(events)
    return events, total


def compare_reports(fortran_text: str, cpp_summary: dict[str, Any], label: str) -> tuple[bool, str]:
    f_events, f_total = parse_fortran_dump(fortran_text)
    c_events = cpp_summary.get("events", [])
    c_total = int(cpp_summary.get("events_count", len(c_events)))

    lines = [
        f"=== Fortran CERNLIB vs DELPHI++ C++20 Cross-Check ({label}) ===",
        f"  Fortran FZIN total events : {f_total}",
        f"  C++20 FZ parser events    : {c_total} (phys_words={cpp_summary.get('physical_record_words')})",
    ]

    ok = True
    if f_total <= 0:
        lines.append("  FAIL: Fortran FZIN reported 0 events.")
        ok = False
    if c_total <= 0:
        lines.append("  FAIL: DELPHI++ C++20 reported 0 events.")
        ok = False
    if f_total != c_total:
        lines.append(f"  FAIL: Event count mismatch (Fortran={f_total} vs C++20={c_total}).")
        ok = False

    n_cmp = min(len(f_events), len(c_events))
    for i in range(n_cmp):
        fe = f_events[i]
        ce = c_events[i]
        lines.append(
            f"  [#{i+1}] F77(run={fe.run}, evt={fe.event}, top={fe.top}, banks={fe.banks}) "
            f"<-> C++(run={ce.get('run')}, evt={ce.get('event')}, banks={ce.get('zebra_banks')}, "
            f"blk={ce.get('raw_blocklets')}, dst={ce.get('dst_particles')}, E_vis={ce.get('visible_energy_gev')})"
        )
        if fe.run > 0 and int(ce.get("run", 0)) > 0 and fe.run != int(ce["run"]):
            lines.append(f"  FAIL: Run mismatch at #{i+1}: Fortran={fe.run} vs C++={ce.get('run')}")
            ok = False

    return ok, "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fortran", type=Path, required=True, help="Fortran dump .txt")
    parser.add_argument("--cpp", type=Path, required=True, help="delphipp_reco --summary-json file")
    parser.add_argument("--label", default="", help="Display label")
    args = parser.parse_args()

    ftext = args.fortran.read_text(encoding="utf-8")
    csum = json.loads(args.cpp.read_text(encoding="utf-8"))
    label = args.label or args.cpp.name
    ok, report = compare_reports(ftext, csum, label)
    print(report)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
