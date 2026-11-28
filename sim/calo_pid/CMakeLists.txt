# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_sim_calo_pid STATIC
  src/calo_pid_sim.cpp
)
add_library(delphipp::sim::calo_pid ALIAS delphipp_sim_calo_pid)

target_include_directories(delphipp_sim_calo_pid
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_sim_calo_pid
  PUBLIC
    delphipp::core
    delphipp::io
)

target_compile_features(delphipp_sim_calo_pid PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_sim_calo_pid PRIVATE /W4)
else()
  target_compile_options(delphipp_sim_calo_pid PRIVATE -Wall -Wextra -Wpedantic)
endif()
