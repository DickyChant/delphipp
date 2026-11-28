# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_reco_calibration STATIC
  src/vd_alignment.cpp
  src/global_alignment.cpp
  src/calo_muon_fixers.cpp
)
add_library(delphipp::reco::calibration ALIAS delphipp_reco_calibration)

target_include_directories(delphipp_reco_calibration
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_reco_calibration
  PUBLIC
    delphipp::core
    delphipp::db
)

target_compile_features(delphipp_reco_calibration PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_reco_calibration PRIVATE /W4)
else()
  target_compile_options(delphipp_reco_calibration PRIVATE -Wall -Wextra -Wpedantic)
endif()
