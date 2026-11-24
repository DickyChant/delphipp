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
