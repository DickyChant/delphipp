# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(delphipp_pipeline STATIC
  src/event_pipeline.cpp
)
add_library(delphipp::pipeline ALIAS delphipp_pipeline)

target_include_directories(delphipp_pipeline
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(delphipp_pipeline
  PUBLIC
    delphipp::core
    delphipp::db
    delphipp::io
    delphipp::dst
    delphipp::reco::hof
    delphipp::reco::mus
    delphipp::reco::muf
    delphipp::reco::mub
    delphipp::reco::tof
    delphipp::reco::emf
    delphipp::reco::had
    delphipp::reco::hpc
    delphipp::reco::stc
    delphipp::reco::sat
    delphipp::reco::vd
    delphipp::reco::tpc
    delphipp::reco::id
    delphipp::reco::od
    delphipp::reco::fwd
    delphipp::reco::rich
    delphipp::reco::tracking
)

target_compile_features(delphipp_pipeline PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_pipeline PRIVATE /W4)
else()
  target_compile_options(delphipp_pipeline PRIVATE -Wall -Wextra -Wpedantic)
endif()
