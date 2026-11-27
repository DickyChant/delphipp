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
| Condition DB & KAPACK | `cargo600.car`, `kapack600.car` | ~30k | `db/` (`ConditionDatabase`, `KapackReader`) | Done | Yes (`ABST0` + `--ddb-kapack`/`load_kapack`) | Binary B-tree & ASCII parity |
| ZEBRA FZ Exchange I/O | `zebra` (`FZ`) | — | `io/` (`FzExchange*`, `fz_zebra_parser`, `raw_blocklet`) | Done | Yes | E2E `FZIN` cross-check |

### Subdetector Reconstruction Modules (All 18 Detectors + RICH)

| Detector | Fortran Package | Non-comment Lines | DELPHI++ Header | Core C++20 Port | Wired in `pipeline/` | Fortran Parity |
|---|---|---:|---|---|---|---|
| **HOF** (Forward Hodoscope) | `hofana.car` | 471 | `reco/hof/hof.hpp` | Done (`unpack_td`) | Yes (`RHOF`/`EGRA`->`TDR`) | Pending full real-data `bank_diff` |
| **MUS** (Surround Muon Chambers) | `musana.car` | 1,233 | `reco/mus/mus.hpp` | Done (`unpack_raw_blocklet`) | Yes (`EDRA`/`MUS`->`TER`) | Pending chamber hit `bank_diff` |
| **MUF** (Forward Muon Chambers) | `mufana.car` | 1,315 | `reco/muf/muf.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RMUF`->`TER`) | Pending chamber hit `bank_diff` |
| **MUB** (Barrel Muon Chambers) | `mubana.car` | 2,783 | `reco/mub/mub.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RMUB`->`TER`) | Pending chamber hit `bank_diff` |
| **TOF** (Time of Flight) | `tofana.car` | 4,259 | `reco/tof/tof.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RTOF`->`TDR`) | Pending timing `bank_diff` |
| **EMF** (Forward EM Calorimeter) | `emfana.car` | 5,181 | `reco/emf/emf.hpp` | Done (`unpack_raw_blocklet`) | Yes (`REMF`->`TER`) | Pending shower `bank_diff` |
| **HAD** (Hadron Calorimeter) | `hadana.car` | 5,878 | `reco/had/had.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RHAC`->`TER`) | Pending tower `bank_diff` |
| **STW** (ID Straw Tube Detector) | `stwana.car` | 6,561 | `reco/stw/stw.hpp` | Done (`unpack_raw_blocklet`) | Yes (`STW`/`STWR`->`TER` + seeding) | Pending real-data straw `bank_diff` |
| **LUM** / **SAT** (VSAT / SAT) | `lumana.car`, `satana.car` | 12,764 | `reco/sat/sat.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RSAT`/`RLUM`->`TER`) | Pending luminosity `bank_diff` |
| **FCA** (Forward Chamber A) | `fcaana.car` | 8,880 | `reco/fwd/fwd_chambers.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RFCA`->`TER`) | Pending wire-plane `bank_diff` |
| **STC** (Small Angle Tile Calo) | `stcana.car` | 11,582 | `reco/stc/stc.hpp` | Done (`unpack_raw_blocklet`) | Yes (`STC`->`TER` + Bhabha) | Pending cluster `bank_diff` |
| **VD** (Vertex Detector) | `vdana.car` | 11,720 | `reco/vd/vd.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RVD `->`TER` + seeding) | Pending alignment DB parity |
| **TPC** (Time Projection Chamber) | `tpcana.car` | 13,368 | `reco/tpc/tpc.hpp` | Done (`94D ABST0`, `dE/dx`, `unpack_raw_blocklet`) | Yes (`RTPC`->`TER` + seeding) | Pending pad/wire `bank_diff` |
| **OD** (Outer Detector) | `odana.car` | 16,834 | `reco/od/od.hpp` | Done (`unpack_raw_blocklet`) | Yes (`ROD `->`TER` + seeding) | Pending plank `bank_diff` |
| **ID** (Inner Detector) | `idana.car` | 17,494 | `reco/id/id.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RID `->`TER` + seeding) | Pending jet-cell `bank_diff` |
| **FCB** (Forward Chamber B) | `fcbana.car` | 18,866 | `reco/fwd/fwd_chambers.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RFCB`->`TER`) | Pending wire-plane `bank_diff` |
| **HPC** (Barrel EM Calorimeter) | `hpcana.car` | 23,187 | `reco/hpc/hpc.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RHPC`->`TER` + Ecal) | Pending 9-layer `bank_diff` |
| **RICH** (Barrel & Forward RICH) | `richall.car` | ~30k | `reco/rich/rich.hpp` | Done (`unpack_raw_blocklet`) | Yes (`RRIB`/`RRIF`->`TER` + PID) | Pending ring likelihood parity |

### Global Tracking, Particle ID, DST, and Simulation

| Subsystem | Fortran Source | Non-comment Lines | DELPHI++ Module | Core C++20 Port | Wired in `pipeline/` | Fortran Parity |
|---|---|---|---|---|---|---|
| Track Ambiguity Solver | `fxsolv.c`, `fxambi.car` | ~5k | `reco/tracking/ambiguity_solver.hpp` | Done | Yes (raw + DST) | Deterministic node budget |
| Kalman Track Refit | `tkfb.car`, `tkff.car`, `tkrefitp.car` | ~18k | `reco/tracking/kalman_fitter.hpp` | Done | Yes (pure-raw seeding -> refit) | Pending multi-surface refit parity |
| Primary & Secondary (`V0`) Vertexing | `pxfvtx.car`, `recv022.car`, `rec2nd.car` | ~12k | `reco/tracking/vertex_fitter.hpp`, `v0_finder.hpp` | Done | Yes (PV in DST; V0 library) | Pending real-data vertex `bank_diff` |
| Energy Flow, Durham Jets & Event Shapes | `eflowxx.car`, `vecsub72.car` | ~8k | `dst/energy_flow.hpp` | Done (`cluster_durham_jets`, thrust, `y_32`) | Library ready | Pending jet `bank_diff` |
| MC Truth Association & `DSTFIX`/`ECORR` | `evasso.car`, `dstfix31.car`, `ecorrxx.car` | ~15k | `dst/event_assoc.hpp`, `dst/dst_fix.hpp` | Done (`associate_particles`, `apply_dst_fixes`) | Library ready | Pending DST fix parity |
| Short/Long DST & PID | `pxdst33.car`, `phdstxx.car`, `shortdst.car`, `muflag.car`, `elepha32.car`, `aabtagxx.car` | ~130k | `dst/` (`short_dst`, `dst_writer`) | Done (`CSV`, `NDJSON`, `HepMC3`) | Yes (raw + `PA`/`ST` -> DST) | E2E `BANK_STAT` + `RAW_STAT` gated |
| Detector Simulation (Phase 4) | `delsim36.car`, `*sim.car` | ~150k | `sim/` *(future phase)* | Not started (out of reco scope) | Container E2E (`runsim`) | E2E Pythia8+DELSIM live |

### What Is Still Not Converted (Among the 164 `.car` Files)
- **Monte Carlo Detector Simulation (`delsim36.car`, `*sim.car`, `geant315.car`, `geanh315.car`)**: Intentional Phase 4 scope — E2E currently runs upstream `runsim` inside `jingyucms/delphi-sim:latest` to generate `.fadsim`/`.sdst` test inputs cleanly in pure C++20 CI.
- **CERNLIB / Legacy MC Generators (`kernlib.car`, `hbook.car`, `jetset*`, `pythi561.car`, `koralz40.car`, `eurodec.car`, `photos.car`, `fatmen.car`, `kupack.car`, `khpack600.car`)**: Replaced natively by C++20 STL, Pythia8, CSV/JSON/HepMC3, and XRootD/EOS filesystem access.
- **Standalone Offline Calibration, Alignment & Strip/Run Fixers (`vdalign.car`, `vdcalyd.car`, `vdbeam.car`, `tpcfix.car`, `fcafix.car`, `fcbfix.car`, `fcbt0.car`, `idbfix.car`, `odbfix.car`, `odtiming.car`, `hacfix36.car`, `mucfix.car`, `ricfix.car`, `moddst.car`, `*cra*.car`)**: Patchy cradles and run-by-run offline calibration production scripts (superseded by `ConditionDatabase::load_kapack` + `apply_dst_fixes`).
- **Full Per-Channel Hardware FASTBUS/LTD Bitfield Tables Across All 18 Subdetectors**: Every subdetector has its working C++20 raw blocklet unpacker and algorithm wired in `pipeline/`; remaining numerical parity tuning replaces generic simulated word masks with exact hardware crate/slot/channel bitfields calibrated against frozen `simana_v94d` DELANA dumps on runs `47014` and `54738`.

### Remaining Engineering Gates to Full 94D Validation
1. **Per-detector hardware FASTBUS/LTD bitfield verification**: verify each subdetector's raw word bitfields against real run `47014`/`54738` DELANA `TDR`/`TER` dumps via `tools/bank_diff.py`.
2. **Offline calibration/alignment table parity**: feed full live `FOR050–056` binary KAPACK geometry (`*DBAS`) and survey constants into every subdetector reconstructor.

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
