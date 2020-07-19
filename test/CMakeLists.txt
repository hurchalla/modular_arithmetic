# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

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
               modular_arithmetic/test_modular_addition.cpp
               modular_arithmetic/test_modular_subtraction.cpp
               modular_arithmetic/test_modular_addition_with_subtraction.cpp
               modular_arithmetic/test_modular_multiplication.cpp
               modular_arithmetic/test_modular_multiplicative_inverse.cpp
               modular_arithmetic/test_modular_pow.cpp
               test_ma.cpp)

EnableMaxWarnings(test_hurchalla_modular_arithmetic)


set_target_properties(test_hurchalla_modular_arithmetic
                      PROPERTIES FOLDER "Tests")
target_link_libraries(test_hurchalla_modular_arithmetic
                      hurchalla_modular_arithmetic
                      gtest_main)
#add_test(test_hurchalla_modular_arithmetic  test_hurchalla_modular_arithmetic)
gtest_discover_tests(test_hurchalla_modular_arithmetic)
