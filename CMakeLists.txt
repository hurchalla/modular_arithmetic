# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET hurchalla_modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.12)
project(hurchalla_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX)

# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

option(TEST_HURCHALLA_LIBS
       "Build tests for all Hurchalla library projects."
       OFF)

# if this is the top level CMakeLists.txt, add testing options, and enable
# testing when testing options have been set to ON.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(TEST_HURCHALLA_TOP_LEVEL
           "Build tests for the top level Hurchalla library project."
           OFF)
    if(TEST_HURCHALLA_TOP_LEVEL OR TEST_HURCHALLA_LIBS)
        enable_testing()
        #include(CTest)
    endif()
elseif(TEST_HURCHALLA_LIBS)
    # If TEST_HURCHALLA_LIBS is set to ON, enable_testing() should have been
    # called either directly or indirectly by the top level project. (Note that
    # if a project calls include(CTest), the included CTest.cmake defines a
    # BUILT_TESTING option and calls enable_testing if BUILD_TESTING is ON.)
    if (NOT CMAKE_TESTING_ENABLED)
        message(FATAL_ERROR "Fatal error: TEST_HURCHALLA_LIBS was set, but enable_testing() was never called")
    endif()
endif()



if(NOT DEFINED CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()


set(public_header_files
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/ma.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/modular_addition.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/modular_multiplication.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/modular_multiplicative_inverse.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/modular_pow.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/modular_subtraction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/traits/type_traits.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/MontgomeryForm.h
    )

set(detail_header_files
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/detail/impl_modular_multiplicative_inverse.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/detail/impl_modular_pow.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/make_safe_unsigned_integer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontgomeryValue.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/monty_common.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontyFullRange.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontySqrtRange.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/sized_uint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h
    )

set(platform_dependent_header_files
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h
    ${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include/hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_multiplication.h
    ${CMAKE_CURRENT_SOURCE_DIR}/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/platform_specific/impl_unsigned_multiply_to_hilo_product.h
    )


#TODO remove  ma.h and .cpp
#TODO I removed a bunch of PBC_WRAP_STDLIB_ASSERT stuff.  make sure contract asserts are wrapped or not, correctly.
#TODO make the asm conditionally included depending on architecture
#TODO make add_library type STATIC or INTERFACE depending on arch


add_library(hurchalla_modular_arithmetic STATIC)

target_sources(hurchalla_modular_arithmetic PRIVATE
            modular_arithmetic/src/ma.cpp
#            modular_arithmetic/src/platform_specific_MSVC_x86_64/modular_multiply_uint64--x64_microsoft.asm
            )

target_include_directories(hurchalla_modular_arithmetic
            PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic/include>)

target_include_directories(hurchalla_modular_arithmetic SYSTEM
            PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)
#use this instead?
#target_include_directories(hurchalla_modular_arithmetic
#            PUBLIC $<INSTALL_INTERFACE:include>)


# Check for required C++ features only if the C++ compiler ID is non-empty.
#target_compile_features(hurchalla_modular_arithmetic INTERFACE 
#          $<$<NOT:$<CXX_COMPILER_ID:>>:
#               cxx_static_assert
#           >)



add_subdirectory(../programming_by_contract
            ${CMAKE_CURRENT_BINARY_DIR}/programming_by_contract)


#add_subdirectory(../common ${CMAKE_CURRENT_BINARY_DIR}/common)


target_link_libraries(hurchalla_modular_arithmetic
                      PRIVATE programming_by_contract)
#target_link_libraries(hurchalla_modular_arithmetic PUBLIC common)
# if you have a public dependency - something your implementation depends upon,
# and that is also visible through your public headers, use PUBLIC above instead
# of PRIVATE




# ***Tests***

# if this is the top level CMakeLists.txt
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    if(TEST_HURCHALLA_TOP_LEVEL OR TEST_HURCHALLA_LIBS)
        add_subdirectory(test)
    endif()
elseif(TEST_HURCHALLA_LIBS)
    add_subdirectory(test)
endif()

