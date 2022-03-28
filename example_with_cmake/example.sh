#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


# This example is meant to show how to use the modular arithmetic library within
# a CMake project.  If you haven't already done so, you should follow the steps
# in the README.md for "How to use the library" | "With CMake"

# Note that using CMake with -DCMAKE_BUILD_TYPE=Release will ensure that the
# standard macro NDEBUG (see <cassert>) is defined, which is essential for best
# performance.

mkdir -p tmp
cmake -S. -B./tmp -DCMAKE_BUILD_TYPE=Release
cmake --build ./tmp --config Release
echo
echo Running example...
echo

./tmp/modular_arithmetic_example
