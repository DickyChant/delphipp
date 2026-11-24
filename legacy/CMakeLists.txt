# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Legacy Fortran reference & strangler-pattern build targets.
# Enabled with -DDELPHIPP_ENABLE_LEGACY_FORTRAN=ON on Linux with gfortran,
# nypatchy/ypatchy, and CERNLIB installed.

if(NOT DELPHIPP_ENABLE_LEGACY_FORTRAN)
  message(STATUS "DELPHI++ legacy Fortran build disabled (enable with -DDELPHIPP_ENABLE_LEGACY_FORTRAN=ON)")
  return()
endif()

find_program(DELPHIPP_PATCHY NAMES nypatchy ypatchy)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

if(NOT EXISTS "${DELPHIPP_SIMANA_DIR}/src/car")
  message(FATAL_ERROR "DELPHIPP_SIMANA_DIR does not contain src/car: ${DELPHIPP_SIMANA_DIR}")
endif()

# Standard gfortran flags for DELPHI / CERNLIB Fortran 77
set(DELPHIPP_LEGACY_FFLAGS
  -std=legacy
  -fno-automatic
  -finit-local-zero
  -fno-range-check
  -fallow-argument-mismatch
  -fd-lines-as-comments
)

# Helper to register a PATCHY cradle expansion + static library target.
# Usage: delphipp_add_legacy_cradle(NAME delana43 CRADLE delana43.cra)
function(delphipp_add_legacy_cradle)
  cmake_parse_arguments(ARG "" "NAME;CRADLE" "EXTRA_FLAGS" ${ARGN})
  set(out_dir "${CMAKE_CURRENT_BINARY_DIR}/expanded/${ARG_NAME}")
  set(out_f "${out_dir}/${ARG_NAME}.f")

  if(DELPHIPP_PATCHY)
    add_custom_command(
      OUTPUT "${out_f}"
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_SOURCE_DIR}/tools/expand_cradle.py"
              --car-dir "${DELPHIPP_SIMANA_DIR}/src/car"
              --cradle "${ARG_CRADLE}"
              --out "${out_f}"
              --patchy "${DELPHIPP_PATCHY}"
      DEPENDS "${DELPHIPP_SIMANA_DIR}/src/car/${ARG_CRADLE}"
              "${CMAKE_SOURCE_DIR}/tools/expand_cradle.py"
      COMMENT "Expanding PATCHY cradle ${ARG_CRADLE} -> ${ARG_NAME}.f"
      VERBATIM
    )
    add_library("legacy_${ARG_NAME}" STATIC "${out_f}")
    target_compile_options("legacy_${ARG_NAME}" PRIVATE ${DELPHIPP_LEGACY_FFLAGS} ${ARG_EXTRA_FLAGS})
    target_link_libraries("legacy_${ARG_NAME}" PUBLIC delphipp::core)
  else()
    message(WARNING "nypatchy/ypatchy not found; skipping legacy target legacy_${ARG_NAME}")
  endif()
endfunction()

# Core Phase 0/1 reference cradles
delphipp_add_legacy_cradle(NAME tanagra322 CRADLE tanagra322.cra)
delphipp_add_legacy_cradle(NAME ux26       CRADLE ux26.cra)
delphipp_add_legacy_cradle(NAME delana43   CRADLE delana43.cra)
delphipp_add_legacy_cradle(NAME delanagen  CRADLE delanagen.cra)
delphipp_add_legacy_cradle(NAME delanamod  CRADLE delanamod.cra)
