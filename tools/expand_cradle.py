#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
"""Expand a DELPHI SIMANA PATCHY cradle (.cra) into Fortran (.f) using nypatchy/ypatchy.

Rewrites legacy directory macros such as (ANADIR), (TANDIR), (UTYDIR), (PHYDIR),
and (DELPHI_PAM) to the local simana_v94d/src/car path and switches platform
flags to LINUX/UNIX before invoking nypatchy.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

DIR_MACROS = (
    "(ANADIR)",
    "(TANDIR)",
    "(UTYDIR)",
    "(PHYDIR)",
    "(DELPHI_PAM)",
    "(SIMDIR)",
)


def rewrite_cradle(text: str, car_dir: Path) -> str:
    car_str = str(car_dir.resolve())
    out_lines: list[str] = ["+USE,LINUX,UNIX."]
    for line in text.splitlines():
        for macro in DIR_MACROS:
            line = line.replace(macro, car_str)
        # Disable VMS/IBM-only +USE lines when targeting Linux
        if re.match(r"^\+USE\s*,\s*\*?VAX", line, re.IGNORECASE):
            line = "+SELF. " + line
        out_lines.append(line)
    return "\n".join(out_lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Expand a DELPHI PATCHY cradle")
    parser.add_argument("--car-dir", type=Path, required=True)
    parser.add_argument("--cradle", type=str, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--patchy", type=str, default="nypatchy")
    args = parser.parse_args()

    cradle_path = args.car_dir / args.cradle
    if not cradle_path.is_file():
        raise SystemExit(f"Cradle not found: {cradle_path}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    rewritten = rewrite_cradle(
        cradle_path.read_text(encoding="latin-1", errors="replace"),
        args.car_dir,
    )

    with tempfile.TemporaryDirectory(prefix="delphipp_patchy_") as tmp:
        tmp_dir = Path(tmp)
        cra_file = tmp_dir / "cradle.cra"
        cra_file.write_text(rewritten, encoding="latin-1")
        out_f = tmp_dir / "out.f"
        log_file = tmp_dir / "patchy.log"

        cmd = [args.patchy, "-", str(out_f), str(cra_file), ".go"]
        proc = subprocess.run(
            cmd,
            cwd=tmp_dir,
            capture_output=True,
            text=True,
            check=False,
        )
        log_file.write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")
        if proc.returncode != 0 or not out_f.is_file():
            raise SystemExit(
                f"{args.patchy} failed (rc={proc.returncode}):\n"
                f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
            )
        shutil.copy2(out_f, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
