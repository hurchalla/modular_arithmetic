# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET hurchalla_basic_modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.13)


# We need to detect if we're using MSVC for x86_64, prior to MSVC2019, since
# these old MSVC versions need a separate asm file for modular multiplication.
# (_MSC_VER < 1920 indicates Visual Studio 2017 or lower)
if((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND (MSVC_VERSION LESS 1920) AND
          (CMAKE_SYSTEM_PROCESSOR MATCHES "x86|X86|amd64|AMD64|EM64T") AND
          (CMAKE_SIZEOF_VOID_P EQUAL 8))
    set(HURCHALLA_MA_USING_OLD_MSVC_X64 TRUE)
else()
    set(HURCHALLA_MA_USING_OLD_MSVC_X64 FALSE)
endif()


if(HURCHALLA_MA_USING_OLD_MSVC_X64)
    project(hurchalla_basic_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX ASM)
else()
    project(hurchalla_basic_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX)
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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_addition.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_multiplication.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_multiplicative_inverse.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_pow.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/modular_subtraction.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/traits/extensible_make_signed.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/traits/extensible_make_unsigned.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/impl_modular_multiplicative_inverse.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/impl_modular_pow.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_multiplication.h>
    )
if(HURCHALLA_MA_USING_OLD_MSVC_X64)
    target_sources(hurchalla_basic_modular_arithmetic PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/platform_specific_MSVC_x86_64/modular_multiply_uint64--x64_microsoft.asm>
        )
endif()


target_include_directories(hurchalla_basic_modular_arithmetic INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)



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
target_link_libraries(hurchalla_basic_modular_arithmetic
                      INTERFACE programming_by_contract)
