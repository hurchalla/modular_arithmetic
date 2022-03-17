# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET test_hurchalla_modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.14)


include(FetchGoogleTest.cmake)
include(EnableMaxWarnings.cmake)


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


# needed for gtest_discover_tests()
include(GoogleTest)


add_executable(test_hurchalla_modular_arithmetic
               modular_arithmetic/test_absolute_value_difference.cpp
               modular_arithmetic/test_modular_addition.cpp
               modular_arithmetic/test_modular_addition_with_subtraction.cpp
               modular_arithmetic/test_modular_multiplication.cpp
               modular_arithmetic/test_modular_multiplicative_inverse.cpp
               modular_arithmetic/test_modular_pow.cpp
               modular_arithmetic/test_modular_subtraction.cpp
               montgomery_arithmetic/low_level_api/test_get_Rsquared_mod_n.cpp
               montgomery_arithmetic/low_level_api/test_inverse_mod_R.cpp
               montgomery_arithmetic/low_level_api/test_REDC.cpp
               montgomery_arithmetic/test_montgomery_pow.cpp
               montgomery_arithmetic/test_MontgomeryForm.cpp
               )

EnableMaxWarnings(test_hurchalla_modular_arithmetic)


set_target_properties(test_hurchalla_modular_arithmetic
                      PROPERTIES FOLDER "Tests")
target_link_libraries(test_hurchalla_modular_arithmetic
                      hurchalla_modular_arithmetic
                      gtest_main)
#add_test(test_hurchalla_modular_arithmetic  test_hurchalla_modular_arithmetic)
gtest_discover_tests(test_hurchalla_modular_arithmetic)
