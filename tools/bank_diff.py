#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
"""ZEBRA / DST bank-level equivalence checker for DELPHI++ validation.

Compares two JSON bank dumps (produced by the Fortran reference and the C++
strangler build) event-by-event and bank-by-bank:
  - Bank structure (bank IDs, link counts, word counts) must match exactly.
  - Integer words must match exactly (unless explicitly masked, e.g. DFIX date/time).
  - Floating-point words must match within configurable absolute/relative tolerances
    (per-bank overrides supported via a JSON policy file).
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class TolerancePolicy:
    default_atol: float = 1e-5
    default_rtol: float = 1e-4
    # Per-bank overrides: bank_id -> {"atol": ..., "rtol": ..., "ignore_words": [1-based idx]}
    bank_overrides: dict[str, dict[str, Any]] = field(default_factory=dict)

    @classmethod
    def from_file(cls, path: Path | None) -> "TolerancePolicy":
        if path is None:
            return cls()
        raw = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            default_atol=float(raw.get("default_atol", 1e-5)),
            default_rtol=float(raw.get("default_rtol", 1e-4)),
            bank_overrides=raw.get("banks", {}),
        )

    def for_bank(self, bank_id: str) -> tuple[float, float, set[int]]:
        ov = self.bank_overrides.get(bank_id.strip(), {})
        atol = float(ov.get("atol", self.default_atol))
        rtol = float(ov.get("rtol", self.default_rtol))
        ignored = {int(x) for x in ov.get("ignore_words", [])}
        return atol, rtol, ignored


def compare_banks(
    ref_events: list[dict[str, Any]],
    cand_events: list[dict[str, Any]],
    policy: TolerancePolicy,
) -> list[str]:
    errors: list[str] = []

    if len(ref_events) != len(cand_events):
        errors.append(
            f"Event count mismatch: ref={len(ref_events)} cand={len(cand_events)}"
        )

    for ev_idx, (ref_ev, cand_ev) in enumerate(zip(ref_events, cand_events)):
        ev_label = f"run={ref_ev.get('run')} evt={ref_ev.get('event', ev_idx)}"
        if ref_ev.get("run") != cand_ev.get("run") or ref_ev.get("event") != cand_ev.get("event"):
            errors.append(
                f"[{ev_idx}] Event header mismatch: ref=({ref_ev.get('run')},{ref_ev.get('event')}) "
                f"cand=({cand_ev.get('run')},{cand_ev.get('event')})"
            )
            continue

        ref_banks = ref_ev.get("banks", [])
        cand_banks = cand_ev.get("banks", [])
        if len(ref_banks) != len(cand_banks):
            errors.append(
                f"[{ev_label}] Bank count mismatch: ref={len(ref_banks)} cand={len(cand_banks)}"
            )
            continue

        for b_idx, (rb, cb) in enumerate(zip(ref_banks, cand_banks)):
            bid = str(rb.get("idh", "")).strip()
            cbid = str(cb.get("idh", "")).strip()
            if bid != cbid or rb.get("idn") != cb.get("idn"):
                errors.append(
                    f"[{ev_label}] Bank #{b_idx} identity mismatch: "
                    f"ref=({bid},{rb.get('idn')}) cand=({cbid},{cb.get('idn')})"
                )
                continue

            atol, rtol, ignored = policy.for_bank(bid)

            # Exact integer word comparison
            ri = rb.get("iwords", {})
            ci = cb.get("iwords", {})
            if set(ri.keys()) != set(ci.keys()):
                errors.append(f"[{ev_label}] Bank {bid} integer word keys differ")
            else:
                for k, rv in ri.items():
                    if int(k) in ignored:
                        continue
                    cv = ci[k]
                    if rv != cv:
                        errors.append(
                            f"[{ev_label}] Bank {bid} int word {k}: ref={rv} cand={cv}"
                        )

            # Tolerant float word comparison
            rf = rb.get("fwords", {})
            cf = cb.get("fwords", {})
            if set(rf.keys()) != set(cf.keys()):
                errors.append(f"[{ev_label}] Bank {bid} float word keys differ")
            else:
                for k, rv in rf.items():
                    if int(k) in ignored:
                        continue
                    cv = cf[k]
                    if not math.isclose(float(rv), float(cv), rel_tol=rtol, abs_tol=atol):
                        errors.append(
                            f"[{ev_label}] Bank {bid} float word {k}: "
                            f"ref={rv} cand={cv} (atol={atol}, rtol={rtol})"
                        )

    return errors


def extract_events(raw: Any) -> list[dict[str, Any]]:
    if isinstance(raw, dict):
        return list(raw.get("events", []))
    if isinstance(raw, list):
        events: list[dict[str, Any]] = []
        for item in raw:
            if isinstance(item, dict) and "events" in item:
                events.extend(item["events"])
            elif isinstance(item, dict):
                events.append(item)
        return events
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare DELPHI ZEBRA/DST bank dumps")
    parser.add_argument("reference", type=Path, help="Golden reference JSON dump")
    parser.add_argument("candidate", type=Path, help="Candidate C++ JSON dump")
    parser.add_argument("--policy", type=Path, default=None, help="Optional tolerance policy JSON")
    args = parser.parse_args()

    ref_events = extract_events(json.loads(args.reference.read_text(encoding="utf-8")))
    cand_events = extract_events(json.loads(args.candidate.read_text(encoding="utf-8")))
    policy = TolerancePolicy.from_file(args.policy)

    errors = compare_banks(ref_events, cand_events, policy)
    if errors:
        print(f"FAIL: {len(errors)} discrepancy(ies) found:")
        for err in errors[:50]:
            print("  -", err)
        if len(errors) > 50:
            print(f"  ... ({len(errors) - 50} more)")
        return 1

    print(f"PASS: {len(ref_events)} event(s) matched within policy.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
