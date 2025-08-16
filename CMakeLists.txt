# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET hurchalla_modular_arithmetic)
    return()
endif()

# later versions are probably fine, but are untested
cmake_minimum_required(VERSION 3.14...4.03)

project(hurchalla_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX)


# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
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


add_library(hurchalla_modular_arithmetic INTERFACE)


add_subdirectory(modular_arithmetic
                 ${CMAKE_CURRENT_BINARY_DIR}/modular_arithmetic)
add_subdirectory(montgomery_arithmetic
                 ${CMAKE_CURRENT_BINARY_DIR}/montgomery_arithmetic)

target_link_libraries(hurchalla_modular_arithmetic
                      INTERFACE hurchalla_basic_modular_arithmetic)
target_link_libraries(hurchalla_modular_arithmetic
                      INTERFACE hurchalla_montgomery_arithmetic)



# TODO:  The following may be overly simple, but works so far to install target
# include directories.  It assumes that the build step from the subdirectories
# montgomery_arithmetic and modular_arithmetic (which have build phase
# target_include_directories commands) provides the information to cmake which
# cmake then uses in the install phase target_include_directories below.  This
# has worked for the basic cmake install tests I've done so far...
# ---------------------
target_include_directories(hurchalla_modular_arithmetic
                           INTERFACE $<INSTALL_INTERFACE:include>)
# TODO: use this instead?
# ----------------
#target_include_directories(hurchalla_modular_arithmetic SYSTEM
#                      INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)



# ***Tests***

# If this is the top level CMakeLists.txt, add testing options, and enable
# testing when testing options have been set to ON.

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(TEST_HURCHALLA_MODULAR_ARITHMETIC
        "Build the tests for the Hurchalla modular arithmetic library project."
        OFF)
    option(HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
        "Include extensive and time consuming tests."
        OFF)
    option(FORCE_TEST_HURCHALLA_CPP11_STANDARD
       "If testing this library, ensure we build googletest and tests using -std=c++11")

    if(TEST_HURCHALLA_MODULAR_ARITHMETIC)
        enable_testing()
        # include(CTest)
        add_subdirectory(test)
    endif()
endif()
