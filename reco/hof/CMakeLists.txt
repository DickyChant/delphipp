# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_reco_hof STATIC
  src/hof.cpp
)
add_library(delphipp::reco::hof ALIAS delphipp_reco_hof)

target_include_directories(delphipp_reco_hof
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_reco_hof PUBLIC delphipp::core)
target_compile_features(delphipp_reco_hof PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_reco_hof PRIVATE /W4)
else()
  target_compile_options(delphipp_reco_hof PRIVATE -Wall -Wextra -Wpedantic)
endif()
