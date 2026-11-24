# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(test_zebra_store test_zebra_store.cpp)
target_link_libraries(test_zebra_store PRIVATE delphipp::core)
add_test(NAME zebra_store COMMAND test_zebra_store)

add_executable(test_hof test_hof.cpp)
target_link_libraries(test_hof PRIVATE delphipp::reco::hof)
add_test(NAME reco_hof COMMAND test_hof)

add_executable(test_condition_db test_condition_db.cpp)
target_link_libraries(test_condition_db PRIVATE delphipp::db)
add_test(NAME db_condition_db COMMAND test_condition_db)

add_executable(test_mus test_mus.cpp)
target_link_libraries(test_mus PRIVATE delphipp::reco::mus)
add_test(NAME reco_mus COMMAND test_mus)

add_executable(test_fz_exchange test_fz_exchange.cpp)
target_link_libraries(test_fz_exchange PRIVATE delphipp::io)
add_test(NAME io_fz_exchange COMMAND test_fz_exchange)

add_executable(test_muf test_muf.cpp)
target_link_libraries(test_muf PRIVATE delphipp::reco::muf)
add_test(NAME reco_muf COMMAND test_muf)

add_executable(test_tanagra_store test_tanagra_store.cpp)
target_link_libraries(test_tanagra_store PRIVATE delphipp::core)
add_test(NAME tanagra_store COMMAND test_tanagra_store)

add_executable(test_mub test_mub.cpp)
target_link_libraries(test_mub PRIVATE delphipp::reco::mub)
add_test(NAME reco_mub COMMAND test_mub)

add_executable(test_field_helix test_field_helix.cpp)
target_link_libraries(test_field_helix PRIVATE delphipp::core)
add_test(NAME core_field_helix COMMAND test_field_helix)

add_executable(test_tof test_tof.cpp)
target_link_libraries(test_tof PRIVATE delphipp::reco::tof)
add_test(NAME reco_tof COMMAND test_tof)

add_executable(test_kapack_reader test_kapack_reader.cpp)
target_link_libraries(test_kapack_reader PRIVATE delphipp::db)
add_test(NAME db_kapack_reader COMMAND test_kapack_reader)

add_executable(test_zebra_bank_io test_zebra_bank_io.cpp)
target_link_libraries(test_zebra_bank_io PRIVATE delphipp::io)
add_test(NAME io_zebra_bank COMMAND test_zebra_bank_io)

add_executable(test_ambiguity_solver test_ambiguity_solver.cpp)
target_link_libraries(test_ambiguity_solver PRIVATE delphipp::reco::tracking)
add_test(NAME reco_ambiguity_solver COMMAND test_ambiguity_solver)

add_executable(test_emf test_emf.cpp)
target_link_libraries(test_emf PRIVATE delphipp::reco::emf)
add_test(NAME reco_emf COMMAND test_emf)

add_executable(test_had test_had.cpp)
target_link_libraries(test_had PRIVATE delphipp::reco::had)
add_test(NAME reco_had COMMAND test_had)

add_executable(test_calorimeters_wave1 test_calorimeters_wave1.cpp)
target_link_libraries(test_calorimeters_wave1 PRIVATE delphipp::reco::hpc delphipp::reco::stc delphipp::reco::sat)
add_test(NAME reco_calorimeters_wave1 COMMAND test_calorimeters_wave1)

add_executable(test_trackers_wave2 test_trackers_wave2.cpp)
target_link_libraries(test_trackers_wave2 PRIVATE delphipp::reco::vd delphipp::reco::tpc delphipp::reco::id delphipp::reco::od delphipp::reco::fwd)
add_test(NAME reco_trackers_wave2 COMMAND test_trackers_wave2)

add_executable(test_rich_vertex_dst test_rich_vertex_dst.cpp)
target_link_libraries(test_rich_vertex_dst PRIVATE delphipp::dst delphipp::reco::rich delphipp::reco::tracking)
add_test(NAME reco_rich_vertex_dst COMMAND test_rich_vertex_dst)

add_executable(test_kalman_and_plumbing test_kalman_and_plumbing.cpp)
target_link_libraries(test_kalman_and_plumbing PRIVATE delphipp::core delphipp::db delphipp::io delphipp::reco::tracking)
add_test(NAME kalman_and_plumbing COMMAND test_kalman_and_plumbing)

add_executable(test_ufield_cargo_dst test_ufield_cargo_dst.cpp)
target_link_libraries(test_ufield_cargo_dst PRIVATE delphipp::core delphipp::db delphipp::dst)
add_test(NAME ufield_cargo_dst COMMAND test_ufield_cargo_dst)

add_executable(test_fz_scanner test_fz_scanner.cpp)
target_link_libraries(test_fz_scanner PRIVATE delphipp::io)
add_test(NAME io_fz_scanner COMMAND test_fz_scanner)

add_executable(test_fz_zebra_parser test_fz_zebra_parser.cpp)
target_link_libraries(test_fz_zebra_parser PRIVATE delphipp::io)
add_test(NAME io_fz_zebra_parser COMMAND test_fz_zebra_parser)
