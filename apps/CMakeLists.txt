# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(delphipp_reco delphipp_reco.cpp)
target_link_libraries(delphipp_reco
  PRIVATE
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
target_compile_features(delphipp_reco PRIVATE cxx_std_20)

if(MSVC)
  target_compile_options(delphipp_reco PRIVATE /W4)
else()
  target_compile_options(delphipp_reco PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(delphipp_fz_inspect delphipp_fz_inspect.cpp)
target_link_libraries(delphipp_fz_inspect PRIVATE delphipp::io)
target_compile_features(delphipp_fz_inspect PRIVATE cxx_std_20)
if(MSVC)
  target_compile_options(delphipp_fz_inspect PRIVATE /W4)
else()
  target_compile_options(delphipp_fz_inspect PRIVATE -Wall -Wextra -Wpedantic)
endif()
