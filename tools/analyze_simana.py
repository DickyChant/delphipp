#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
"""Static analyzer for DELPHI SIMANA PATCHY (.car) source trees.

Extracts:
  - inventory.json : per-file metrics (lines, decks, routines, idioms)
  - commons.json   : COMMON block -> files & routines mapping
  - callgraph.json : routine -> called routines (static CALL edges)
  - banks.json     : ZEBRA/TANAGRA bank creation & CDE (+KEEP) catalogue
  - SUMMARY.md     : human-readable summary & low-coupling pilot candidates
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

RE_PATCH = re.compile(r"^\+PATCH\s*,\s*([^.,\s]+)", re.IGNORECASE)
RE_DECK = re.compile(r"^\+DECK\s*,\s*([^.,\s]+)", re.IGNORECASE)
RE_KEEP = re.compile(r"^\+KEEP\s*,\s*([^.,\s]+)", re.IGNORECASE)
RE_CDE = re.compile(r"^\+(?:CDE|SEQ)\s*,\s*([^.,\s]+)", re.IGNORECASE)
RE_SELF_IF = re.compile(r"^\+SELF\s*,.*IF\s*=", re.IGNORECASE)
RE_ROUTINE = re.compile(
    r"^\s{6}\s*(?:SUBROUTINE|(?:[A-Z0-9*]+\s+)?FUNCTION|PROGRAM|BLOCK\s*DATA)\s+([A-Z0-9_]+)",
    re.IGNORECASE,
)
RE_COMMON = re.compile(r"\bCOMMON\s*/\s*([A-Z0-9_]+)\s*/", re.IGNORECASE)
RE_CALL = re.compile(r"\bCALL\s+([A-Z0-9_]+)", re.IGNORECASE)
RE_MZBOOK = re.compile(
    r"\bCALL\s+(MZBOOK|MZLIFT|MZDROP|MZGARB|MZWIPE|MZPUSH)\b", re.IGNORECASE
)
RE_ZEBRA_ARR = re.compile(r"\b(LQ|IQ|Q)\s*\(", re.IGNORECASE)
RE_GOTO = re.compile(r"\bGO\s*TO\b", re.IGNORECASE)
RE_EQUIV = re.compile(r"\bEQUIVALENCE\b", re.IGNORECASE)
RE_ENTRY = re.compile(r"^\s{6}\s*ENTRY\s+([A-Z0-9_]+)", re.IGNORECASE)


def analyze_directory(car_dir: Path) -> dict:
    car_files = sorted(car_dir.glob("*.car"))
    file_records: list[dict] = []
    commons_map: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"files": set(), "routines": set(), "keeps": set()}
    )
    callgraph: dict[str, set[str]] = defaultdict(set)
    keeps_catalogue: dict[str, dict] = {}
    mz_calls: dict[str, list[dict]] = defaultdict(list)

    for car_path in car_files:
        text = car_path.read_text(encoding="latin-1", errors="replace")
        lines = text.splitlines()

        current_patch = ""
        current_deck = ""
        current_keep = ""
        current_routine = f"<file:{car_path.stem}>"

        rec = {
            "file": car_path.name,
            "stem": car_path.stem.lower(),
            "total_lines": len(lines),
            "code_lines": 0,
            "patchy_directives": 0,
            "patches": 0,
            "decks": 0,
            "keeps": 0,
            "cde_includes": 0,
            "self_if_guards": 0,
            "routines": [],
            "entries": 0,
            "commons": set(),
            "zebra_array_refs": 0,
            "mz_calls": 0,
            "gotos": 0,
            "equivalences": 0,
        }

        for lineno, line in enumerate(lines, 1):
            if not line.strip():
                continue
            first = line[0]
            if first == "+":
                rec["patchy_directives"] += 1
                if m := RE_PATCH.match(line):
                    current_patch = m.group(1).upper()
                    current_keep = ""
                    rec["patches"] += 1
                elif m := RE_DECK.match(line):
                    current_deck = m.group(1).upper()
                    current_keep = ""
                    rec["decks"] += 1
                elif m := RE_KEEP.match(line):
                    current_keep = m.group(1).upper()
                    rec["keeps"] += 1
                    keeps_catalogue.setdefault(
                        current_keep,
                        {
                            "defined_in": car_path.name,
                            "patch": current_patch,
                            "line": lineno,
                            "commons": set(),
                            "used_in_files": set(),
                        },
                    )
                elif m := RE_CDE.match(line):
                    rec["cde_includes"] += 1
                    seq = m.group(1).upper()
                    keeps_catalogue.setdefault(
                        seq,
                        {
                            "defined_in": None,
                            "patch": "",
                            "line": 0,
                            "commons": set(),
                            "used_in_files": set(),
                        },
                    )["used_in_files"].add(car_path.name)
                if RE_SELF_IF.match(line):
                    rec["self_if_guards"] += 1
                continue

            if first in ("C", "c", "*", "!"):
                continue

            rec["code_lines"] += 1

            if m := RE_ROUTINE.match(line):
                current_routine = m.group(1).upper()
                rec["routines"].append(current_routine)
            elif m := RE_ENTRY.match(line):
                rec["entries"] += 1

            for m in RE_COMMON.finditer(line):
                cname = m.group(1).upper()
                rec["commons"].add(cname)
                commons_map[cname]["files"].add(car_path.name)
                if current_keep:
                    commons_map[cname]["keeps"].add(current_keep)
                    keeps_catalogue[current_keep]["commons"].add(cname)
                else:
                    commons_map[cname]["routines"].add(current_routine)

            for m in RE_CALL.finditer(line):
                callee = m.group(1).upper()
                callgraph[current_routine].add(callee)

            if m := RE_MZBOOK.search(line):
                rec["mz_calls"] += 1
                mz_calls[m.group(1).upper()].append(
                    {
                        "file": car_path.name,
                        "routine": current_routine,
                        "line": lineno,
                    }
                )

            rec["zebra_array_refs"] += len(RE_ZEBRA_ARR.findall(line))
            if RE_GOTO.search(line):
                rec["gotos"] += 1
            if RE_EQUIV.search(line):
                rec["equivalences"] += 1

        rec["commons"] = sorted(rec["commons"])
        rec["routine_count"] = len(rec["routines"])
        file_records.append(rec)

    commons_out = {
        k: {
            "file_count": len(v["files"]),
            "files": sorted(v["files"]),
            "routine_count": len(v["routines"]),
            "keeps": sorted(v["keeps"]),
        }
        for k, v in sorted(commons_map.items())
    }

    keeps_out = {
        k: {
            "defined_in": v["defined_in"],
            "patch": v["patch"],
            "line": v["line"],
            "commons": sorted(v["commons"]),
            "used_in_file_count": len(v["used_in_files"]),
        }
        for k, v in sorted(keeps_catalogue.items())
    }

    callgraph_out = {k: sorted(v) for k, v in sorted(callgraph.items())}

    return {
        "inventory": file_records,
        "commons": commons_out,
        "callgraph": callgraph_out,
        "banks": {
            "keeps": keeps_out,
            "mz_call_counts": {k: len(v) for k, v in sorted(mz_calls.items())},
        },
    }


def write_summary(data: dict, out_path: Path) -> None:
    inv = data["inventory"]
    tot_lines = sum(r["total_lines"] for r in inv)
    tot_code = sum(r["code_lines"] for r in inv)
    tot_routines = sum(r["routine_count"] for r in inv)
    tot_commons = len(data["commons"])

    det_ana = [
        r
        for r in inv
        if r["stem"].endswith("ana")
        and r["stem"] not in ("delana43", "dstana10", "hctanaxx")
    ]
    det_ana.sort(key=lambda r: (r["code_lines"], len(r["commons"])))

    top_commons = sorted(
        data["commons"].items(),
        key=lambda kv: kv[1]["file_count"],
        reverse=True,
    )[:20]

    lines = [
        "# SIMANA v94D Static Inventory Summary",
        "",
        f"- **`.car` files:** {len(inv)}",
        f"- **Total lines:** {tot_lines:,} ({tot_code:,} non-comment code lines)",
        f"- **Routines:** {tot_routines:,}",
        f"- **Distinct COMMON blocks:** {tot_commons:,}",
        f"- **Distinct PATCHY `+KEEP` / `+CDE` sequences:** {len(data['banks']['keeps']):,}",
        "",
        "## Phase 2 Pilot Candidates (Subdetector `*ana.car`, smallest first)",
        "",
        "| Package | Code lines | Routines | COMMONs | ZEBRA refs | GOTOs |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in det_ana:
        lines.append(
            f"| `{r['file']}` | {r['code_lines']:,} | {r['routine_count']} | "
            f"{len(r['commons'])} | {r['zebra_array_refs']:,} | {r['gotos']:,} |"
        )

    lines += [
        "",
        "## Most Widely Shared COMMON Blocks (Top 20 by file count)",
        "",
        "| COMMON | Files | Defined in `+KEEP` |",
        "|---|---:|---|",
    ]
    for cname, info in top_commons:
        keeps = ", ".join(info["keeps"][:4]) if info["keeps"] else "*(inline)*"
        lines.append(f"| `/{cname}/` | {info['file_count']} | {keeps} |")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze DELPHI SIMANA .car files")
    parser.add_argument("car_dir", type=Path, help="Path to simana_v94d/src/car")
    parser.add_argument("--out", type=Path, required=True, help="Output directory")
    args = parser.parse_args()

    data = analyze_directory(args.car_dir)
    args.out.mkdir(parents=True, exist_ok=True)

    (args.out / "inventory.json").write_text(
        json.dumps(data["inventory"], indent=2) + "\n", encoding="utf-8"
    )
    (args.out / "commons.json").write_text(
        json.dumps(data["commons"], indent=2) + "\n", encoding="utf-8"
    )
    (args.out / "callgraph.json").write_text(
        json.dumps(data["callgraph"], indent=2) + "\n", encoding="utf-8"
    )
    (args.out / "banks.json").write_text(
        json.dumps(data["banks"], indent=2) + "\n", encoding="utf-8"
    )
    write_summary(data, args.out / "SUMMARY.md")
    print(f"Wrote inventory to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
