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
