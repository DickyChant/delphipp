<!--
SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
SPDX-License-Identifier: GPL-3.0-or-later
-->
# DELPHI++ (`delphipp`)

Pure-C++20 standalone rewrite of the **DELPHI** simulation and reconstruction software (`simana_v94d`, targeting 1994D LEP1 parity), with **zero Fortran or CERNLIB runtime dependencies**.

- **API Documentation (GitHub Pages):** https://dickychant.github.io/delphipp/
- **Full Reconstruction Roadmap:** [docs/ROADMAP.md](docs/ROADMAP.md)
- **Static Fortran Inventory Summary:** [docs/inventory/SUMMARY.md](docs/inventory/SUMMARY.md)

## Architecture Overview

```text
DELSIM / DELANA ZEBRA FZ (.fadsim / .sdst / .fz)
  │
  ▼
io/fz_zebra_parser   (BE/LE Exchange physical blocks, steering/fast blocks, bank materializer)
  │
  ▼
core/ZebraStoreView  +  db/ConditionDatabase (94D .asc + binary FOR050-056 KAPACK B-tree)
  │                     core/MagneticField   (real 41-barrel + 39-endcap UFIDAT solenoid map)
  ▼
pipeline/EventReconstructionPipeline
  ├── Raw blocklets & subdetector unpacking (18 detectors -> TanagraStore TDR/TER)
  ├── Global tracking: ambiguity solver (fxsolv) -> Kalman refit (Highland MS) -> TKR
  └── Primary vertex fit + Particle ID (TPC dE/dx, RICH, Calorimetry, Muon chambers)
  │
  ▼
dst/ShortDST  ──►  CSV  |  NDJSON  |  HepMC3 Asciiv3  |  bank_diff.py JSON
                                                            ▲
                                                            │  (CI E2E cross-check)
tools/fortran_fz_dump.f (CERNLIB FZIN reference in DELPHI container)
```

## Module Replacement Roadmap (`simana_v94d` -> C++20)

Static analysis of `simana_v94d` (`164` `.car` files, `810,903` non-comment lines, `10,418` routines, `2,735` COMMON blocks) defines the exact replacement surface below. Status is tracked across three explicit gates:

1. **Core C++20 Port**: Pure-C++20 algorithms, data structures, and unit tests implemented.
2. **Wired in `pipeline/`**: Connected to `EventReconstructionPipeline` and `delphipp_reco` on real ZEBRA FZ files.
3. **Fortran Parity**: Exact per-event/per-bank equivalence validated against CERNLIB Fortran references on runs `47014` and `54738` (and self-contained 10-event Pythia8/DELSIM E2E).

### Core Infrastructure, I/O, and Condition Database

| Subsystem | Fortran Source | Non-comment Lines | DELPHI++ Module | Core C++20 Port | Wired in `pipeline/` | Fortran Parity |
|---|---|---:|---|---|---|---|
| ZEBRA / UX Store | `ux26.car` | ~10k | `core/` (`ZebraStoreView`, `OwnedZebraStore`) | Done | Yes |FZ structural match |
| TANAGRA Bank Store | `tanagra322.car` | ~10k | `core/` (`TanagraStore`: `TDR`/`TER`/`TKR`) | Done | Yes | Schema ready |
| Geometry Transforms | `ddapp600.car` | ~13k | `core/` (`dtrans`) | Done | Yes | Unit tested |
| Magnetic Field | `ufield.car` | ~1.5k | `core/` (`MagneticField`, `UFIDAT` map) | Done | Yes | Grid verified |
| Helix Extrapolation | `exx.car` | ~8k | `core/` (`PerigeeTrack` cylinder/z-plane) | Done | Yes | Unit tested |
| Condition DB & KAPACK | `cargo600.car`, `kapack600.car` | ~30k | `db/` (`ConditionDatabase`, `KapackReader`) | Done | Yes (`ABST0`) | Binary DB decoding |
| ZEBRA FZ Exchange I/O | `zebra` (`FZ`) | — | `io/` (`FzExchange*`, `fz_zebra_parser`, `raw_blocklet`) | Done | Yes | E2E `FZIN` cross-check |

### Subdetector Reconstruction Modules (All 18 Detectors + RICH)

| Detector | Fortran Package | Non-comment Lines | DELPHI++ Header | Core C++20 Port | Wired in `pipeline/` | Fortran Parity |
|---|---|---:|---|---|---|---|
| **HOF** (Forward Hodoscope) | `hofana.car` | 471 | `reco/hof/hof.hpp` | Done | Yes (`EGRA`/`EQRA`) | Pending full real-data `bank_diff` |
| **MUS** (Surround Muon Chambers) | `musana.car` | 1,233 | `reco/mus/mus.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **MUF** (Forward Muon Chambers) | `mufana.car` | 1,315 | `reco/muf/muf.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **MUB** (Barrel Muon Chambers) | `mubana.car` | 2,783 | `reco/mub/mub.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **TOF** (Time of Flight) | `tofana.car` | 4,259 | `reco/tof/tof.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **SAT** (Small Angle Tagger) | `satana.car` | 5,100 | `reco/sat/sat.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **EMF** (Forward EM Calorimeter) | `emfana.car` | 5,181 | `reco/emf/emf.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **HAD** (Hadron Calorimeter) | `hadana.car` | 5,878 | `reco/had/had.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **STW** / **LUM** (Straw / Lumi) | `stwana.car`, `lumana.car` | 14,225 | `reco/stc/stc.hpp`, `reco/sat/sat.hpp` | Partial (combined in `stc`/`sat`) | Blocklet decode | Pending standalone split & `bank_diff` |
| **FCA** (Forward Chamber A) | `fcaana.car` | 8,880 | `reco/fwd/fwd_chambers.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **STC** (Small Angle Tile Calo) | `stcana.car` | 11,582 | `reco/stc/stc.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **VD** (Vertex Detector) | `vdana.car` | 11,720 | `reco/vd/vd.hpp` | Done | Blocklet decode | Pending cluster/alignment DB parity |
| **TPC** (Time Projection Chamber) | `tpcana.car` | 13,368 | `reco/tpc/tpc.hpp` | Done (`94D ABST0`, `dE/dx`) | Blocklet + `ABST0` | Pending raw pad/wire unpacker parity |
| **OD** (Outer Detector) | `odana.car` | 16,834 | `reco/od/od.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **ID** (Inner Detector) | `idana.car` | 17,494 | `reco/id/id.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **FCB** (Forward Chamber B) | `fcbana.car` | 18,866 | `reco/fwd/fwd_chambers.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **HPC** (Barrel EM Calorimeter) | `hpcana.car` | 23,187 | `reco/hpc/hpc.hpp` | Done | Blocklet decode | Pending raw unpacker & `bank_diff` |
| **RICH** (Barrel & Forward RICH) | `richall.car` | ~30k | `reco/rich/rich.hpp` | Done | Blocklet + PID | Pending photon ring likelihood parity |

### Global Tracking, Particle ID, DST, and Simulation

| Subsystem | Fortran Source | Non-comment Lines | DELPHI++ Module | Core C++20 Port | Wired in `pipeline/` | Fortran Parity |
|---|---|---:|---|---|---|---|
| Track Ambiguity Solver | `fxsolv.c`, `fxambi.car` | ~5k | `reco/tracking/ambiguity_solver.hpp` | Done | Yes | Deterministic node budget |
| Kalman Track Refit | `tkfb.car`, `tkff.car`, `tkrefitp.car` | ~18k | `reco/tracking/kalman_fitter.hpp` | Done | Yes | Pending multi-surface refit parity |
| Primary Vertex Fit | `pxfvtx.car` | ~6k | `reco/tracking/vertex_fitter.hpp` | Done | Yes | Pending real-data vertex `bank_diff` |
| Short/Long DST & PID | `pxdst33.car`, `phdstxx.car`, `shortdst.car`, `muflag.car`, `elepha32.car`, `aabtagxx.car` | ~130k | `dst/` (`short_dst`, `dst_writer`) | Done (`CSV`, `NDJSON`, `HepMC3`) | Yes (`PA`/`PV` + PID) | Real `.sdst` decoded in E2E |
| Detector Simulation (Phase 4) | `delsim36.car`, `sx*.car` | ~150k | `sim/` *(future phase)* | Planned | Container E2E (`runsim`) | E2E Pythia8+DELSIM live |

### Remaining Engineering Gates to Full v1.0 94D Parity
1. **Per-detector FASTBUS/LTD raw-to-TD unpackers**: expand `pipeline/` beyond `HOF` (`EGRA`) + generic blocklet headers so every subdetector decodes raw FASTBUS/LTD words into full physical hits (`VD`, `ID`, `TPC`, `OD`, `FCA`/`FCB`, `HPC`, `EMF`, `HAD`, `STC`, `SAT`, `TOF`, `MUB`/`MUF`/`MUS`, `RICH`).
2. **Real binary condition DB tables**: wire `db::KapackReader` directly to DELPHI `FOR050–056` binary files for geometry/calibration constants (`*DBAS`).
3. **Full pattern recognition -> Kalman -> Vertex -> Short-DST**: drive `ShortDst` purely from unpacked raw `.fadsim` hits (without needing pre-built `.sdst` `PA`/`PV` banks), matching Fortran `TDR`/`TER`/`TKR`/`PA`/`PV` outputs via `tools/bank_diff.py` on runs `47014` and `54738`.
4. **Full bank-level Fortran-vs-C++ tolerance gating**: extend `tools/compare_fortran_cpp.py` from event/header/bank-count checks to full `TDR`/`TER`/`PA`/`PV` tolerance checks against Fortran DELANA outputs.

## Repository Layout

```text
delphipp/
├── core/        # C++20 ZEBRA store, TanagraStore, MagneticField, Helix, dtrans
├── db/          # Condition database (ASCII .asc + binary FOR050-056 KAPACK B-tree)
├── io/          # ZEBRA FZ Exchange reader/writer, scanner, bank-tree parser, raw blocklets
├── reco/        # All 18 subdetectors + global tracking (ambiguity, Kalman, vertex)
├── dst/         # Short-DST builder + CSV / NDJSON / HepMC3 Asciiv3 writers
├── pipeline/    # Unified EventReconstructionPipeline orchestrating all layers
├── apps/        # delphipp_reco & delphipp_fz_inspect CLI binaries
├── tools/       # analyze_simana.py, bank_diff.py, fortran_fz_dump.f, compare_fortran_cpp.py
├── tests/       # C++20 unit and integration tests (CTest)
├── docs/        # ROADMAP.md, Doxyfile (doxygen-awesome-css), static inventory
└── docker/      # Reproducible Linux toolchain
```

## Quick Start (Linux / Container)

```bash
# Build C++20 libraries, apps, & run full CTest suite
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build
ctest --test-dir build --output-on-failure

# Run standalone reconstruction on real DELSIM/DELANA .fadsim or .sdst FZ files
./build/apps/delphipp_reco \
  --input simana.sdst \
  --verbose \
  --out reco_banks.json \
  --summary-json reco_summary.json \
  --csv-out dst_particles.csv \
  --hepmc-out events.hepmc3
```

## License

Licensed under **GPL-3.0-or-later** (see [LICENSES/GPL-3.0-or-later.txt](LICENSES/GPL-3.0-or-later.txt)), consistent with `simana_v94d`.
