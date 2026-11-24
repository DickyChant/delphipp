<!--
SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
SPDX-License-Identifier: GPL-3.0-or-later
-->
# DELPHI++ (`delphipp`)

Incremental C++20 rewrite of the DELPHI simulation and reconstruction software (`simana_v94d`, LEP1 1994D processing).

## Strategy

DELPHI++ follows the **strangler pattern**: instead of a clean-slate rewrite, we build a mixed Fortran/C++ executable where legacy Fortran packages are replaced by C++20 implementations one package at a time behind compatible ABIs, validated after every swap against frozen Fortran reference DSTs (runs `47014` and `54738`).

## Repository Layout

```text
delphipp/
├── core/        # C++20 ZEBRA-compatible store, UX/TANAGRA bridge, event loop
├── db/          # Condition database snapshot reader (FOR050–FOR056)
├── io/          # ZEBRA FZ I/O (raw, full DST, short DST) and future ROOT/RNTuple
├── reco/        # Subdetector reconstruction & global tracking (Phase 2)
├── sim/         # DELSIM detector simulation (Phase 4)
├── dst/         # PXDST, PHDST, short/long DST production (Phase 3)
├── legacy/      # PATCHY expansion & CMake targets for frozen Fortran reference
├── tools/       # Static analyzer, PATCHY expander, bank-diff validation harness
├── tests/       # C++ unit tests and regression test drivers
├── docs/        # Generated codebase inventory, bank catalogues, architecture notes
└── docker/      # Reproducible Linux toolchain (C++20 + gfortran + CERNLIB)
```

## Quick Start (Linux / Container)

```bash
# Build C++20 core & unit tests
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build
ctest --test-dir build --output-on-failure

# Build with legacy Fortran reference (requires gfortran, CERNLIB, and simana_v94d)
cmake -B build-legacy -G Ninja \
  -DDELPHIPP_ENABLE_LEGACY_FORTRAN=ON \
  -DDELPHIPP_SIMANA_DIR=../simana_v94d
cmake --build build-legacy
```

## Static Analysis & Inventory

```bash
# Python 3 (Linux/macOS/Windows)
python3 tools/analyze_simana.py ../simana_v94d/src/car --out docs/inventory

# PowerShell (Windows, zero dependencies)
pwsh tools/analyze_simana.ps1 -CarDir ..\simana_v94d\src\car -OutDir docs\inventory
```

## License

Licensed under **GPL-3.0-or-later** (see [LICENSES/GPL-3.0-or-later.txt](LICENSES/GPL-3.0-or-later.txt)), consistent with `simana_v94d`.
