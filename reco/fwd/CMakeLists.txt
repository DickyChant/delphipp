# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_reco_fwd STATIC src/fwd_chambers.cpp)
add_library(delphipp::reco::fwd ALIAS delphipp_reco_fwd)
target_include_directories(delphipp_reco_fwd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(delphipp_reco_fwd PUBLIC delphipp::core)
target_compile_features(delphipp_reco_fwd PUBLIC cxx_std_20)
if(MSVC)
  target_compile_options(delphipp_reco_fwd PRIVATE /W4)
else()
  target_compile_options(delphipp_reco_fwd PRIVATE -Wall -Wextra -Wpedantic)
endif()
