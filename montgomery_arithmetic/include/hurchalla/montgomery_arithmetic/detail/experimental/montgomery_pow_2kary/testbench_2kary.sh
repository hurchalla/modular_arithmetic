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


#cppcompiler=g++
#cppcompiler=clang++
cppcompiler=$1

#optimization_level=O2
#optimization_level=O3
optimization_level=$2

#define_mont_type=-DDEF_MONT_TYPE=MontgomeryQuarter<U>
define_mont_type=-DDEF_MONT_TYPE=$3
define_uint_type=-DDEF_UINT_TYPE=$4


# argument $8 (if present), should be -DHURCHALLA_ALLOW_INLINE_ASM_ALL
define_use_asm=$8


cpp_standard=c++17


# You need to clone the util, factoring, and modular_arithmetic repos
# from https://github.com/hurchalla


# SET repo_directory TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT
# REPOSITORIES.  (or otherwise ensure the compiler /I flags correctly specify
# the needed hurchalla include directories)
repo_directory=/Users/jeffreyhurchalla/Desktop



if [[ $cppcompiler == "g++" ]]; then
  error_limit=-fmax-errors=3
  warn_nrvo=-Wnrvo
else
  error_limit=-ferror-limit=3
fi



# argument $8 (if present), should be -DHURCHALLA_ALLOW_INLINE_ASM_ALL


# to debug you can compile with the below also
# -DHURCHALLA_CLOCKWORK_ENABLE_ASSERTS  -DHURCHALLA_UTIL_ENABLE_ASSERTS

# we could also use  -g  to get debug symbols (for lldb/gdb, and objdump)


$cppcompiler  \
        $error_limit   -$optimization_level \
        $define_mont_type  $define_uint_type  $define_use_asm \
        -Wall -Wextra -Wpedantic $warn_nrvo  \
        -std=$cpp_standard \
        -I${repo_directory}/modular_arithmetic/modular_arithmetic/include \
        -I${repo_directory}/modular_arithmetic/montgomery_arithmetic/include \
        -I${repo_directory}/util/include \
        -c testbench_montgomery_pow_2kary.cpp

$cppcompiler  -$optimization_level  -std=$cpp_standard  -o testbench_montgomery_pow_2kary  testbench_montgomery_pow_2kary.o -lm



exit_on_failure

echo "compilation finished, now executing:"


# argument $5 (if present), is the randomization seed for std::mt19937_64
# argument $6 (if present), is max_modulus_bits_reduce
# argument $7 (if present), is exponent_bits_reduce

./testbench_montgomery_pow_2kary $5 $6 $7

# To give you an example of invoking this script at the command line:
#   ./testbench.sh clang++ O3 MontgomeryFull __uint128_t 191 8 50 -DHURCHALLA_ALLOW_INLINE_ASM_ALL


