# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET hurchalla_montgomery_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.13)

project(hurchalla_montgomery_arithmetic VERSION 1.0.0 LANGUAGES CXX)


# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


add_library(hurchalla_montgomery_arithmetic INTERFACE)


# Check for required C++ features only if the C++ compiler ID is non-empty.
#target_compile_features(hurchalla_montgomery_arithmetic INTERFACE 
#          $<$<NOT:$<CXX_COMPILER_ID:>>:
#               cxx_static_assert
#           >)


target_sources(hurchalla_montgomery_arithmetic INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/MontgomeryForm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontgomeryValue.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/monty_common.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyFullRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontySqrtRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/sized_uint.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/impl_unsigned_multiply_to_hilo_product.h>
    )

target_include_directories(hurchalla_montgomery_arithmetic
            INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)


add_subdirectory(../modular_arithmetic
            ${CMAKE_CURRENT_BINARY_DIR}/modular_arithmetic)
target_link_libraries(hurchalla_montgomery_arithmetic
                      INTERFACE hurchalla_basic_modular_arithmetic)


# TODO: Currently we assume the programming_by_contract dependency is in a
# folder at the relative location specified below.  This is a hack, since users
# usually won't know this is a requirement, and it obligates the user to perform
# a manual step (even if only once) to get the dependency from git and place it
# in the correct location.  The dependency should be completely invisible to the
# user, and this cmake script should "just work" without any extra steps.
# Among the possible solutions might be git submodules, or FetchContent.  Git
# submodules often seems to be warned against, making me more optimistic about
# FetchContent or perhaps some other solution.  Note: the tests for this project
# currently use FetchContent to get the Google Tests dependency.
# ---------------------------
add_subdirectory(../../programming_by_contract
            ${CMAKE_CURRENT_BINARY_DIR}/programming_by_contract)
target_link_libraries(hurchalla_montgomery_arithmetic
                      INTERFACE programming_by_contract)
