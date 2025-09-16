# Copyright (c) 2020-2025 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET hurchalla_montgomery_arithmetic)
    return()
endif()

# later versions are probably fine, but are untested
cmake_minimum_required(VERSION 3.14...4.03)

project(hurchalla_montgomery_arithmetic VERSION 1.0.0 LANGUAGES CXX)


# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


add_library(hurchalla_montgomery_arithmetic INTERFACE)


target_sources(hurchalla_montgomery_arithmetic INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/MontgomeryForm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/montgomery_form_aliases.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/BaseMontgomeryValue.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/ImplMontgomeryForm.contents>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/ImplMontgomeryForm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/impl_montgomery_two_pow.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontgomeryFormExtensions.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyFullRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyTags.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/montgomery_pow_2kary.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/MontyFullRangeMasked.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryForm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryWrapper.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/ConcreteMontgomeryForm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/montgomery_pow.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/montgomery_two_pow.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/two_times_restricted.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/quarterrange_get_canonical.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/detail/platform_specific/halfrange_get_canonical.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/get_Rsquared_mod_n.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/REDC.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/detail/impl_inverse_mod_R.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_array_get_Rsquared_mod_n.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_get_Rsquared_mod_n.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/ImplRedc.h>
    )


install(DIRECTORY
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    DESTINATION include)


target_include_directories(hurchalla_montgomery_arithmetic
            INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)


add_subdirectory(../modular_arithmetic
            ${CMAKE_CURRENT_BINARY_DIR}/modular_arithmetic)
target_link_libraries(hurchalla_montgomery_arithmetic
                      INTERFACE hurchalla_basic_modular_arithmetic)


include(FetchContent)
FetchContent_Declare(
    hurchalla_util
    GIT_REPOSITORY https://github.com/hurchalla/util.git
    GIT_TAG        22181fce5ae399896921a8f0538a88a9dbe41e66
)
FetchContent_MakeAvailable(hurchalla_util)

target_link_libraries(hurchalla_montgomery_arithmetic
                      INTERFACE hurchalla_util)
