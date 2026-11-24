# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_reco_tracking STATIC
  src/ambiguity_solver.cpp
  src/vertex_fitter.cpp
)
add_library(delphipp::reco::tracking ALIAS delphipp_reco_tracking)

target_include_directories(delphipp_reco_tracking
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_reco_tracking PUBLIC delphipp::core)
target_compile_features(delphipp_reco_tracking PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_reco_tracking PRIVATE /W4)
else()
  target_compile_options(delphipp_reco_tracking PRIVATE -Wall -Wextra -Wpedantic)
endif()
