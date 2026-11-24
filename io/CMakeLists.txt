# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_io STATIC
  src/fz_exchange.cpp
  src/fz_scanner.cpp
  src/raw_blocklet.cpp
  src/zebra_bank_io.cpp
)
add_library(delphipp::io ALIAS delphipp_io)

target_include_directories(delphipp_io
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_io PUBLIC delphipp::core)
target_compile_features(delphipp_io PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_io PRIVATE /W4)
else()
  target_compile_options(delphipp_io PRIVATE -Wall -Wextra -Wpedantic)
endif()
