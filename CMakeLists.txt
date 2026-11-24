# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.22)

project(delphipp
  VERSION 0.1.0
  DESCRIPTION "Incremental C++20 rewrite of DELPHI SIMANA v94D"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(DELPHIPP_BUILD_TESTS "Build C++ unit and integration tests" ON)
option(DELPHIPP_ENABLE_LEGACY_FORTRAN "Enable mixed Fortran/C++ strangler build against simana_v94d" OFF)
set(DELPHIPP_SIMANA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../simana_v94d" CACHE PATH "Path to simana_v94d checkout")

if(DELPHIPP_ENABLE_LEGACY_FORTRAN)
  enable_language(Fortran)
  enable_language(C)
endif()

add_subdirectory(core)
add_subdirectory(db)
add_subdirectory(io)
add_subdirectory(reco)
add_subdirectory(dst)
add_subdirectory(apps)
add_subdirectory(legacy)

if(DELPHIPP_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
