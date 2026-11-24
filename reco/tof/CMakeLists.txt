# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_reco_tof STATIC
  src/tof.cpp
)
add_library(delphipp::reco::tof ALIAS delphipp_reco_tof)

target_include_directories(delphipp_reco_tof
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_reco_tof PUBLIC delphipp::core)
target_compile_features(delphipp_reco_tof PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_reco_tof PRIVATE /W4)
else()
  target_compile_options(delphipp_reco_tof PRIVATE -Wall -Wextra -Wpedantic)
endif()
