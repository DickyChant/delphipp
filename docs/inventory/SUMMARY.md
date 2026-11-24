# SIMANA v94D Static Inventory Summary

- **`.car` files:** 164
- **Total lines:** 1,265,354 (810,903 non-comment code lines)
- **Routines:** 10,418
- **Distinct COMMON blocks:** 2,735
- **Distinct PATCHY `+KEEP` / `+CDE` sequences:** 3,273

## Phase 2 Pilot Candidates (Subdetector `*ana.car`, smallest first)

| Package | Code lines | Routines | COMMONs | ZEBRA refs | GOTOs |
|---|---:|---:|---:|---:|---:|
| `hofana.car` | 471 | 10 | 4 | 26 | 40 |
| `musana.car` | 1,233 | 17 | 4 | 6 | 31 |
| `mufana.car` | 1,315 | 18 | 6 | 24 | 76 |
| `mubana.car` | 2,783 | 34 | 5 | 56 | 125 |
| `tofana.car` | 4,259 | 40 | 27 | 106 | 87 |
| `satana.car` | 5,100 | 35 | 32 | 160 | 173 |
| `emfana.car` | 5,181 | 77 | 13 | 903 | 236 |
| `hadana.car` | 5,878 | 84 | 25 | 252 | 409 |
| `stwana.car` | 6,561 | 75 | 23 | 31 | 280 |
| `lumana.car` | 7,664 | 48 | 44 | 287 | 182 |
| `fcaana.car` | 8,880 | 74 | 14 | 313 | 512 |
| `stcana.car` | 11,582 | 172 | 41 | 387 | 366 |
| `vdana.car` | 11,720 | 115 | 20 | 75 | 241 |
| `tpcana.car` | 13,368 | 126 | 53 | 1,339 | 428 |
| `odana.car` | 16,834 | 182 | 51 | 56 | 627 |
| `idana.car` | 17,494 | 196 | 33 | 812 | 594 |
| `fcbana.car` | 18,866 | 189 | 105 | 156 | 653 |
| `hpcana.car` | 23,187 | 247 | 101 | 2,393 | 840 |

## Most Widely Shared COMMON Blocks (Top 20 by file count)

| COMMON | Files | Defined in `+KEEP` |
|---|---:|---|
| `/PAWC/` | 23 | EXHMEM, HBOOK_MEMORY, HCBOOK, PAWC |
| `/QUEST/` | 22 | FATBANK, HB_AXI, Q, QC_COMMON |
| `/LUJETS/` | 13 | LUJEDA, LUJETS, LUNDCOM |
| `/SLATE/` | 10 | ALL, SLATE |
| `/ZUNIT/` | 10 | UXZLUN, ZUNIT |
| `/CFREAD/` | 9 | CFREAD, EXFFRD, SXFFRD, UXCFF |
| `/LUDAT1/` | 9 | LUDAT1, LUNDCOM |
| `/LUDAT3/` | 8 | LUDAT3, LUJEDA, LUNDCOM |
| `/DDAPCM/` | 8 | DDAMEM, DDSPACE, EFCALBUN, EXCDDAPP |
| `/LUDAT2/` | 7 | LUDAT2, LUJEDA, LUNDCOM |
| `/ZTELLC/` | 7 | *(inline)* |
| `/HEPEVT/` | 6 | HEPEVT |
| `/LUDAT4/` | 6 | LUDAT4, LUNDCOM |
| `/FZSTAT/` | 5 | FZSTAT, PHLUN, UXFZSTAT |
| `/PXCONS/` | 5 | PXCONS |
| `/LUDATR/` | 5 | LUDATR, LUDATR1, LUDATR2 |
| `/PXCHDR/` | 5 | PXCHDR |
| `/ZUNITZ/` | 5 | ZUNIT |
| `/PYSUBS/` | 4 | PYSUBS |
| `/RKZFCO/` | 4 | RKZFCO |
