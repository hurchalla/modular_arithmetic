#!/bin/bash

# Copyright (c) 2025 Jeffrey Hurchalla.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


exit_on_failure () {
  if [ $? -ne 0 ]; then
    exit 1
  fi
}


cppcompiler=clang++
cpp_standard="-std=c++17"


# You need to clone the util, factoring, and modular_arithmetic repos
# from https://github.com/hurchalla


# SET repo_directory TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT
# REPOSITORIES.  (or otherwise ensure the compiler /I flags correctly specify
# the needed hurchalla include directories)
repo_directory=/home/jeff/repos


$cppcompiler  \
        -O3  -DNDEBUG \
        $cpp_standard \
        -I${repo_directory}/modular_arithmetic/modular_arithmetic/include \
        -I${repo_directory}/modular_arithmetic/montgomery_arithmetic/include \
        -I${repo_directory}/util/include \
        -c example_montgomery_two_pow.cpp

$cppcompiler  -O3  -std="c++17"  -o example_montgomery_two_pow  example_montgomery_two_pow.o -lm



exit_on_failure

echo "compilation finished, now executing:"


./example_montgomery_two_pow
