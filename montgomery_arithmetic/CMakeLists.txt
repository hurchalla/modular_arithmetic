# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET hurchalla_montgomery_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.14)

project(hurchalla_montgomery_arithmetic VERSION 1.0.0 LANGUAGES CXX)


# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


add_library(hurchalla_montgomery_arithmetic INTERFACE)


target_sources(hurchalla_montgomery_arithmetic INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/MontgomeryForm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/optimization_tag_structs.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/inverse_mod_r.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/monty_common.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyFullRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontySixthRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/RedcSmallR.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/sized_uint.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/impl_unsigned_multiply_to_hilo_product.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/MontHelper.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/RedcLargeR.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/MontySqrtRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/platform_specific/montadd_sqrt_range.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/platform_specific/montsub_sqrt_range.h>
    )

target_include_directories(hurchalla_montgomery_arithmetic
            INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)


add_subdirectory(../modular_arithmetic
            ${CMAKE_CURRENT_BINARY_DIR}/modular_arithmetic)
target_link_libraries(hurchalla_montgomery_arithmetic
                      INTERFACE hurchalla_basic_modular_arithmetic)


include(FetchContent)
FetchContent_Declare(
    programming_by_contract
    GIT_REPOSITORY https://github.com/hurchalla/programming_by_contract.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(programming_by_contract)

target_link_libraries(hurchalla_montgomery_arithmetic
                      INTERFACE hurchalla_programming_by_contract)
