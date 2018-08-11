# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET test_modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.12)
#project(test_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX)


include(FetchGoogleTest.cmake)


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


# needed for gtest_discover_tests()
include(GoogleTest)


add_executable(test_modular_arithmetic  test_ma.cpp)
set_target_properties(test_modular_arithmetic PROPERTIES FOLDER "Tests")
target_link_libraries(test_modular_arithmetic
                      modular_arithmetic
                      gtest_main)
#add_test(test_modular_arithmetic  test_modular_arithmetic)
gtest_discover_tests(test_modular_arithmetic)
