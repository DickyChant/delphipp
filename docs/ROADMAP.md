<!--
SPDX-FileCopyrightText: 2026 Sitian Qian
SPDX-License-Identifier: GPL-3.0-or-later
-->
# DELPHI++ Full Reconstruction Roadmap

**North-Star Goal:** Standalone, pure-C++20 DELPHI event reconstruction (**DELANA + Full/Short DST**, targeting 1994D reprocessing parity) with zero Fortran or CERNLIB runtime dependencies, validated against frozen `simana_v94d` reference outputs on runs `47014` and `54738`.

## Milestone Gates

| Milestone | Scope | Status |
|---|---|---|
| **v0.1 — Phase 0 Scaffold** | CMake/C++20, CI (GCC + Clang + Docker), static inventory (`164` `.car` files), `bank_diff.py` harness | Done |
| **v0.2 — Vertical Slice** | `ZebraStoreView`, `TanagraStore`, `DTRANS`, `MagneticField`, `PerigeeTrack` helix, ASCII `ConditionDatabase`, `FzExchange` framing, `HOF`/`MUS`/`MUF`/`MUB`/`TOF` core reco, `delphipp_reco` CLI + CI `bank_diff` gate | Done |
| **v0.3 — Real Binary I/O & Tracking Core** | Binary `FOR050–056` KAPACK reader (`db/`), ZEBRA FZ bank materializer (`io/`), global track ambiguity solver (`reco/tracking/`) | In progress |
| **v0.4 — Real-Data Outer/Muon Validation** | Complete raw unpackers (`*DECO`/`*CALF`/`*UNPB`) & DB loaders (`*DBAS`) for `HOF`/`MUS`/`MUF`/`MUB`/`TOF`; match Fortran `TDR`/`TER` banks on real 94D raw events | Planned |
| **v0.5 — Calorimeters & Luminosity** | `SAT`, `EMF`, `HAC`, `STW`, `LUM`, `STC`, `HPC` | Planned |
| **v0.6 — Tracking Detectors & Global Fit** | `VD`, `ID`, `TPC`, `OD`, `FCA`, `FCB`, Kalman refit (`tkfb`/`tkff`/`exx`), vertexing (`vertex`/`fitvdz`) | Planned |
| **v0.7 — RICH & Particle ID** | `richall`, `dedx*`, `muflag`, `elepha32`, `aabtagxx` | Planned |
| **v1.0 — Full 94D DELANA + DST Parity** | `pxdst33`, `phdstxx`, `dstana10`, `shortdst`, `longdst`, `*fix` decks, full validation on runs `47014` & `54738` | Planned |

## Subsystem Progress Matrix

| Layer | Subsystem | Fortran Source | Non-comment Lines | Pure C++20 Status |
|---|---|---|---:|---|
| Core | ZEBRA & UX store | `ux26.car` | ~10k | `ZebraStoreView` + `OwnedZebraStore` + UX COMMONs |
| Core | TANAGRA bank store | `tanagra322.car` | ~10k | `TanagraStore` (`TDR`/`TER`/`TKR` + `bank_diff` JSON) |
| Core | Geometry transforms | `ddapp600.car` | ~13k | `dtrans()` |
| Core | Magnetic field | `ufield.car` | ~1.5k | `MagneticField` solenoid + SAT quad (grid tables pending) |
| Core | Helix extrapolation | `exx.car` | ~8k | `PerigeeTrack` analytical cylinder & z-plane propagation |
| DB | Condition DB & KAPACK | `cargo600.car`, `kapack600.car` | ~30k | ASCII `.asc` store + `KapackReader` binary block/header reader |
| I/O | ZEBRA FZ Exchange | `zebra` | — | `FzExchangeWriter`/`Reader` + bank materializer |
| Reco | HOF | `hofana.car` | 471 | `reco/hof` (`EGTD`, `EGTE`, `EGNBUN`) |
| Reco | MUS | `musana.car` | 1,233 | `reco/mus` (`EDAMSH`, `EDCLSR`, `EDCBCL`, `EDBLTE`) |
| Reco | MUF | `mufana.car` | 1,315 | `reco/muf` (`EQTSLW`, `EQDLDF`, `EQWPOS`, `EQDRIF`, `EQSPAF`) |
| Reco | MUB | `mubana.car` | 2,783 | `reco/mub` (`EMCHAN`, `EMSLPB`, `EMCHID`, `EMLTIM`, `EMDFTT`) |
| Reco | TOF | `tofana.car` | 4,259 | `reco/tof` (`EFSPA`, `EFNBUN`) |
| Reco | Tracking ambiguity | `fxsolv.c`, `fxambi.car` | ~5k | `reco/tracking` ambiguity solver |
| Reco | Calorimeters (7) | `sat`/`emf`/`had`/`stw`/`lum`/`stc`/`hpc` | ~66k | Planned (v0.5) |
| Reco | Trackers (6) + RICH | `vd`/`tpc`/`od`/`id`/`fca`/`fcb`/`richall` | ~111k | Planned (v0.6–v0.7) |
| DST | PXDST / PHDST / Short DST | `pxdst33`/`phdstxx`/`shortdst`/... | ~130k | Planned (v1.0) |
