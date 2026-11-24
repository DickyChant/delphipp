# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(test_zebra_store test_zebra_store.cpp)
target_link_libraries(test_zebra_store PRIVATE delphipp::core)
add_test(NAME zebra_store COMMAND test_zebra_store)

add_executable(test_hof test_hof.cpp)
target_link_libraries(test_hof PRIVATE delphipp::reco::hof)
add_test(NAME reco_hof COMMAND test_hof)
