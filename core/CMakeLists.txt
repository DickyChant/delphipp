# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_core STATIC
  src/tanagra_store.cpp
  src/zebra_store.cpp
)
add_library(delphipp::core ALIAS delphipp_core)

target_include_directories(delphipp_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(delphipp_core PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_core PRIVATE /W4)
else()
  target_compile_options(delphipp_core PRIVATE -Wall -Wextra -Wpedantic)
endif()
