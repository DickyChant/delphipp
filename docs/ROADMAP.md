<!--
SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
SPDX-License-Identifier: GPL-3.0-or-later
-->
# DELPHI++ Full Reconstruction Roadmap

**North-Star Goal:** Standalone, pure-C++20 DELPHI event reconstruction (**DELANA + Full/Short DST**, targeting 1994D reprocessing parity) with zero Fortran or CERNLIB runtime dependencies, validated against frozen `simana_v94d` reference outputs on runs `47014` and `54738` and live 10-event DELSIM E2E runs.

See [README.md](../README.md) for the live module-by-module replacement matrix (published at https://dickychant.github.io/delphipp/).

## Milestone Gates

| Milestone | Scope | Status |
|---|---|---|
| **v0.1 — Phase 0 Scaffold** | CMake/C++20, CI (GCC + Clang + Docker), static inventory (`164` `.car` files), `bank_diff.py` harness | Done |
| **v0.2 — Vertical Slice** | `ZebraStoreView`, `TanagraStore`, `dtrans`, `MagneticField`, `PerigeeTrack` helix, ASCII `ConditionDatabase`, `FzExchange` framing, `HOF`/`MUS`/`MUF`/`MUB`/`TOF` core reco, `delphipp_reco` CLI + CI `bank_diff` gate | Done |
| **v0.3 — Real Binary I/O & Tracking Core** | Binary `FOR050–056` KAPACK reader (`db/`), ZEBRA FZ Exchange scanner & bank-tree parser (`io/`), global track ambiguity solver (`reco/tracking/`) | Done |
| **v0.4 — Calorimeters & Luminosity** | Pure-C++20 algorithms for `EMF`, `HAD`, `HPC`, `STC`, `SAT`/`LUM` | Done (core C++20) |
| **v0.5 — Tracking Detectors, RICH & Vertex** | Pure-C++20 algorithms for `VD`, `ID`, `TPC` (94D `ABST0` & truncated `dE/dx`), `OD`, `FCA`/`FCB`, `RICH`, Kalman refit, primary vertex fitter | Done (core C++20) |
| **v0.6 — Unified Pipeline, Short DST & E2E** | `pipeline/` (`EventReconstructionPipeline`), `dst/` (`ShortDST`, `CSV`/`NDJSON`/`HepMC3`), `fortran_fz_dump.f` + `compare_fortran_cpp.py`, Doxygen Awesome + GitHub Pages, 10-event Pythia8/DELSIM E2E | Done |
| **v0.7 — Full Raw Blocklet Unpackers & Binary DDB** | Decode raw FASTBUS/LTD blocklets into full `TDR`/`TER` hits across all 18 detectors and load live constants from `FOR050–056` via `KapackReader` | Next |
| **v1.0 — Full 94D DELANA + DST Parity** | Drive `ShortDST` solely from raw `.fadsim` hits (full tracking + PID without `.sdst` pre-built `PA`/`PV`), and match frozen Fortran `simana_v94d` banks on runs `47014` & `54738` within `bank_diff.py` tolerances | In progress |

## Static Codebase Scale (`simana_v94d`)

- **`.car` files:** `164` (`810,903` non-comment Fortran lines, `10,418` routines, `2,735` COMMON blocks)
- Full per-package inventory & dependency graph: [docs/inventory/SUMMARY.md](inventory/SUMMARY.md)
