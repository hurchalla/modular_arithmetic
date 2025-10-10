# Copyright (c) 2020-2025 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET hurchalla_basic_modular_arithmetic)
    return()
endif()

# later versions are probably fine, but are untested
cmake_minimum_required(VERSION 3.14...4.03)

project(hurchalla_basic_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX)

# We need to detect if we're using MSVC for x86_64, prior to MSVC2019, since
# these old MSVC versions need a separate asm file for modular multiplication.
# (_MSC_VER < 1920 indicates Visual Studio 2017 or lower)
if((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND (MSVC_VERSION LESS 1920) AND
          (CMAKE_SYSTEM_PROCESSOR MATCHES "x86|X86|amd64|AMD64|EM64T") AND
          (CMAKE_SIZEOF_VOID_P EQUAL 8))
    set(HURCHALLA_MA_USING_OLD_MSVC_X64 TRUE)
    enable_language(ASM_MASM)
else()
    set(HURCHALLA_MA_USING_OLD_MSVC_X64 FALSE)
endif()


# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


if(HURCHALLA_MA_USING_OLD_MSVC_X64)
    add_library(hurchalla_basic_modular_arithmetic STATIC)
else()
    add_library(hurchalla_basic_modular_arithmetic INTERFACE)
endif()


target_sources(hurchalla_basic_modular_arithmetic INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/absolute_value_difference.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_addition.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_multiplication.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_multiplicative_inverse.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_pow.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_subtraction.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/impl_modular_multiplicative_inverse.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/impl_modular_pow.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/optimization_tag_structs.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/platform_specific/impl_absolute_value_difference.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_addition.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_multiplication.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_subtraction.h>
    )
if(HURCHALLA_MA_USING_OLD_MSVC_X64)
    target_sources(hurchalla_basic_modular_arithmetic PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/platform_specific_MSVC_x86_64/modular_multiply_uint64--x64_microsoft.asm>
        )
endif()


install(DIRECTORY
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    DESTINATION include)


target_include_directories(hurchalla_basic_modular_arithmetic INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)


include(FetchContent)
FetchContent_Declare(
    hurchalla_util
    GIT_REPOSITORY https://github.com/hurchalla/util.git
    GIT_TAG        5ff4af01d0a769968e56bc668f4c9a9ad1896749
)
FetchContent_MakeAvailable(hurchalla_util)

target_link_libraries(hurchalla_basic_modular_arithmetic
                      INTERFACE hurchalla_util)
