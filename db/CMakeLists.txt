# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_db STATIC
  src/condition_db.cpp
)
add_library(delphipp::db ALIAS delphipp_db)

target_include_directories(delphipp_db
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(delphipp_db PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_db PRIVATE /W4)
else()
  target_compile_options(delphipp_db PRIVATE -Wall -Wextra -Wpedantic)
endif()
