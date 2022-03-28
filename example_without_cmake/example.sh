#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


# This example is intended for the case that you are not using CMake.
# If you haven't already done so, you should follow the steps in the README.md
# for "How to use the library" | "Without CMake"


# ------------------------------------------------------------------------
# You'll need to change the include_path below, and you may need to change
# the cpp_compiler.
# ------------------------------------------------------------------------

# set include_path to the directory where you installed the modular arithmetic
# library, with "/include" appended to the end of that directory name.
include_path=/home/jeff/Desktop/include

# set the compiler to whatever you want.  Below is gcc or clang.
#cpp_compiler=g++
cpp_compiler=clang++


$cpp_compiler -std="c++17" \
        -Wall -Wextra -ferror-limit=3  \
        -O2  -DNDEBUG  \
        -I$include_path \
        -o example  example.cpp

./example  
