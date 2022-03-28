# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


cmake_minimum_required(VERSION 3.14)

project(modular_arithmetic_example VERSION 1.0.0 LANGUAGES CXX)

add_executable(modular_arithmetic_example)

target_sources(modular_arithmetic_example PRIVATE
    example.cpp
    )

# in add_subdirectory below, the first argument must be the path on your system
# to the root of the cloned modular arithmetic repository.  In this example,
# ".." is used because (unless you move it) this CMakeLists.txt file has a path
# that is within the repository, and that path is one directory level below the
# repo's root.
add_subdirectory(..  ${CMAKE_CURRENT_BINARY_DIR}/modular_arithmetic)

target_link_libraries(modular_arithmetic_example
                      hurchalla_modular_arithmetic)



# To build this example, see the file example.sh.
