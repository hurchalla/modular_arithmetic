# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET test_hurchalla_modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.13)


include(FetchGoogleTest.cmake)


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


# needed for gtest_discover_tests()
include(GoogleTest)


add_executable(test_hurchalla_modular_arithmetic
               test_ma.cpp)


# Set compiler warnings.  We want our library headers to compiler cleanly
# with as many warnings as possible, so that we know users will not get
# compiler warnings from our headers.
if(MSVC)
    target_compile_options(test_hurchalla_modular_arithmetic PRIVATE /W4 /WX)
else()
    target_compile_options(test_hurchalla_modular_arithmetic PRIVATE -Werror
          -Wall -Wextra -Wpedantic -pedantic-errors -Wshadow -Wcast-qual
          -Wmissing-include-dirs -Wredundant-decls -Wnon-virtual-dtor
          -Wstrict-overflow=4 -Wconversion -Wsign-conversion -Wfloat-equal
          -Winvalid-pch -Wstack-protector -Wpacked -Wnull-dereference
          -Wuninitialized -Wregister -Wunused -Wunused-parameter -Wformat=2
          -Wdisabled-optimization -Wvla -Wold-style-cast -Wdouble-promotion
          -Winline -Wundef -Wformat-nonliteral -Wmissing-noreturn
          -Wzero-as-null-pointer-constant -Weffc++ -Wctor-dtor-privacy
          -Wmissing-declarations)
    if((NOT (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")) OR
               (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.1))
        target_compile_options(test_hurchalla_modular_arithmetic PRIVATE
                -Wextra-semi)
    endif()

    # additional compiler specific warnings
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(test_hurchalla_modular_arithmetic PRIVATE
                -ferror-limit=3  -Wcast-align -Wmismatched-tags)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(test_hurchalla_modular_arithmetic PRIVATE
                -fmax-errors=3 -Wlogical-op -Wnoexcept
                -Wduplicated-cond -Wduplicated-branches
                -Wstringop-overflow=3 -Wsuggest-final-methods -Wplacement-new=2
                -Wsuggest-final-types -Wsuggest-override -Walloc-zero
                -Wnormalized -Wformat-truncation=2 -Wformat-overflow=2
                -Wformat-signedness -Wstrict-null-sentinel -Wopenmp-simd
                -Wunsafe-loop-optimizations -Wvector-operation-performance)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.1)
            target_compile_options(test_hurchalla_modular_arithmetic PRIVATE
                               -Wcast-align=strict -Wcomma-subscript -Wvolatile)
        else()
            target_compile_options(test_hurchalla_modular_arithmetic PRIVATE
                               -Wcast-align)
        endif()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # TODO: investigate icc warnings
    endif()
endif()


set_target_properties(test_hurchalla_modular_arithmetic
                      PROPERTIES FOLDER "Tests")
target_link_libraries(test_hurchalla_modular_arithmetic
                      hurchalla_modular_arithmetic
                      gtest_main)
#add_test(test_hurchalla_modular_arithmetic  test_hurchalla_modular_arithmetic)
gtest_discover_tests(test_hurchalla_modular_arithmetic)
