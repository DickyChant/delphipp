# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_dst STATIC
  src/dst_fix.cpp
  src/dst_writer.cpp
  src/energy_flow.cpp
  src/event_assoc.cpp
  src/short_dst.cpp
)
add_library(delphipp::dst ALIAS delphipp_dst)
target_include_directories(delphipp_dst PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(delphipp_dst PUBLIC delphipp::core delphipp::reco::rich delphipp::reco::tracking)
target_compile_features(delphipp_dst PUBLIC cxx_std_20)
if(MSVC)
  target_compile_options(delphipp_dst PRIVATE /W4)
else()
  target_compile_options(delphipp_dst PRIVATE -Wall -Wextra -Wpedantic)
endif()
