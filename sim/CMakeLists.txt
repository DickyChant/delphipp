# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_sim STATIC
  src/sim_core.cpp
  src/tracker_sim.cpp
)
add_library(delphipp::sim ALIAS delphipp_sim)

target_include_directories(delphipp_sim
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_sim
  PUBLIC
    delphipp::core
    delphipp::io
    delphipp::reco::vd
    delphipp::reco::id
    delphipp::reco::tpc
    delphipp::reco::od
    delphipp::reco::fwd
    delphipp::reco::stw
)

target_compile_features(delphipp_sim PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_sim PRIVATE /W4)
else()
  target_compile_options(delphipp_sim PRIVATE -Wall -Wextra -Wpedantic)
endif()
