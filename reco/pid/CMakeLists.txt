# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_reco_pid STATIC
  src/electron_id.cpp
  src/muon_id.cpp
  src/hadron_id.cpp
  src/btag.cpp
)
add_library(delphipp::reco::pid ALIAS delphipp_reco_pid)

target_include_directories(delphipp_reco_pid
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_reco_pid
  PUBLIC
    delphipp::core
    delphipp::reco::rich
    delphipp::reco::tpc
    delphipp::reco::tracking
)
target_compile_features(delphipp_reco_pid PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_reco_pid PRIVATE /W4)
else()
  target_compile_options(delphipp_reco_pid PRIVATE -Wall -Wextra -Wpedantic)
endif()
