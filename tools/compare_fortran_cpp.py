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
    r"EVT\s+#\s*(\d+)\s+run=\s*(-?\d+)\s+evt=\s*(-?\d+)\s+top=(.{4}?)\s+nd=\s*(\d+)\s+banks=\s*(\d+)"
)
BANK_STAT_RE = re.compile(
    r"BANK_STAT\s+npa=\s*(\d+)\s+nchg=\s*(\d+)\s+nneu=\s*(\d+)\s+sum_e=\s*([-+0-9.eE]+)\s+sum_p=\s*([-+0-9.eE]+)"
)
RAW_STAT_RE = re.compile(r"RAW_STAT\s+nrawd=\s*(\d+)\s+npmb=\s*(\d+)")
TOTAL_RE = re.compile(r"Total events read by CERNLIB FZIN:\s*(\d+)")
NON_EVENT_BANKS = {"RUN", "CORR", "RUNH"}


@dataclass
class FortranEvent:
    index: int
    run: int
    event: int
    top: str
    nd: int
    banks: int
    npa: int | None = None
    nchg: int | None = None
    nneu: int | None = None
    sum_e: float | None = None
    sum_p: float | None = None
    nrawd: int | None = None
    npmb: int | None = None


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
                    top=m.group(4).strip(),
                    nd=int(m.group(5)),
                    banks=int(m.group(6)),
                )
            )
            continue
        if events:
            mb = BANK_STAT_RE.search(line)
            if mb:
                events[-1].npa = int(mb.group(1))
                events[-1].nchg = int(mb.group(2))
                events[-1].nneu = int(mb.group(3))
                events[-1].sum_e = float(mb.group(4))
                events[-1].sum_p = float(mb.group(5))
                continue
            mr = RAW_STAT_RE.search(line)
            if mr:
                events[-1].nrawd = int(mr.group(1))
                events[-1].npmb = int(mr.group(2))
    mt = TOTAL_RE.search(text)
    total = int(mt.group(1)) if mt else len(events)
    return events, total


def compare_reports(fortran_text: str, cpp_summary: dict[str, Any], label: str) -> tuple[bool, str]:
    f_all, f_raw_total = parse_fortran_dump(fortran_text)
    f_phys = [e for e in f_all if e.top not in NON_EVENT_BANKS]
    f_other = [e for e in f_all if e.top in NON_EVENT_BANKS]
    f_total = len(f_phys) if f_all else f_raw_total

    c_events = cpp_summary.get("events", [])
    c_total = int(cpp_summary.get("events_count", len(c_events)))

    lines = [
        f"=== Fortran CERNLIB vs DELPHI++ C++20 Cross-Check ({label}) ===",
        f"  Fortran FZIN records      : total={f_raw_total}, physics_and_sim={f_total}, preamble={len(f_other)}",
        f"  C++20 FZ parser events    : {c_total} (phys_words={cpp_summary.get('physical_record_words')})",
    ]
    if f_other:
        lines.append("  Fortran preamble banks    : " + ", ".join(f"{e.top}(nd={e.nd})" for e in f_other))

    ok = True
    if f_total <= 0:
        lines.append("  FAIL: Fortran FZIN reported 0 physics/sim events.")
        ok = False
    if c_total <= 0:
        lines.append("  FAIL: DELPHI++ C++20 reported 0 events.")
        ok = False
    if f_total != c_total:
        lines.append(f"  FAIL: Record count mismatch (Fortran={f_total} vs C++20={c_total}).")
        ok = False

    n_cmp = min(len(f_phys), len(c_events))
    for i in range(n_cmp):
        fe = f_phys[i]
        ce = c_events[i]
        c_run = int(ce.get("run", 0))
        c_evt = int(ce.get("event", -1))
        c_top = str(ce.get("top_bank", "")).strip() or "?"
        c_dets = ce.get("detectors", {})
        det_str = ",".join(f"{k}:{v}" for k, v in sorted(c_dets.items()))
        det_suffix = f", dets=[{det_str}]" if det_str else ""
        lines.append(
            f"  [#{i+1}] F77(run={fe.run}, evt={fe.event}, top={fe.top}, banks={fe.banks}) "
            f"<-> C++(run={c_run}, evt={c_evt}, top={c_top}, banks={ce.get('zebra_banks')}, "
            f"blk={ce.get('raw_blocklets')}, dst={ce.get('dst_particles')}, E_vis={ce.get('visible_energy_gev')}{det_suffix})"
        )
        if fe.run != 0 and c_run != 0 and abs(fe.run) != abs(c_run):
            lines.append(f"  FAIL: Run mismatch at #{i+1}: Fortran={fe.run} vs C++={c_run}")
            ok = False
        if fe.event != c_evt:
            lines.append(f"  FAIL: Event mismatch at #{i+1}: Fortran={fe.event} vs C++={c_evt}")
            ok = False
        if fe.npa is not None:
            c_dst = int(ce.get("dst_particles", 0))
            if fe.npa != c_dst:
                lines.append(
                    f"  FAIL: DTOP PA->DST particle count mismatch at #{i+1}: Fortran npa={fe.npa} vs C++ dst={c_dst}"
                )
                ok = False
            if fe.sum_e is not None:
                c_evis = float(ce.get("visible_energy_gev", 0.0))
                tol = max(2.0, 0.05 * fe.sum_e)
                if abs(fe.sum_e - c_evis) > tol:
                    lines.append(
                        f"  FAIL: DTOP visible energy mismatch at #{i+1}: Fortran={fe.sum_e:.3f} vs C++={c_evis:.3f} (tol={tol:.3f})"
                    )
                    ok = False
        if fe.nrawd is not None:
            c_blk = int(ce.get("raw_blocklets", 0))
            if fe.nrawd != c_blk:
                lines.append(
                    f"  FAIL: TOP RAWD blocklet count mismatch at #{i+1}: Fortran nrawd={fe.nrawd} vs C++ blk={c_blk}"
                )
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
