#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


# Description of this script -----------
# This is a working convenience script for invoking the testing builds and then
# running the tests.
# The syntax is 
# ./build_tests [-c<compiler_name>] [-r] [-a] [-u] [-s] [-m<Release|Debug>]
#
# -c allows you to select the compiler, rather than using the default.
# -r specifies to run all tests after the build.  Without -r, no tests will run.
# -a specifies you want to compile the code using typically helpful (how much it
#    helps depends on your compiler) inline asm optimizations, which makes for
#    the fastest binaries but of course has the downsides of inline asm -
#    primarily that inline asm is extremely difficult to properly test.
# -u specifies that you want to compile the code using all available inline asm
#    routines, so that the tests will cover all of them (this is not expected to
#    result in the fastest binaries).
# -s specifies that you wish to compile the code with optimizations using a
#    series of bitwise operations, instead of operations that typically would
#    compile to conditional select or conditional move instructions.  This
#    usually improves performance on architectures that lack those instructions;
#    in certain cases it can also improve performance for architectures that
#    possess those instructions.
# -m allows you to choose between Release and Debug build configuration, rather
#    than using the default.
#
# Currently it supports clang, gcc, and icc but you'll need to customize the
# section under "#Compiler commands" to match the compilers on your system.  The
# script doesn't do anything fancy like auto detection or anything like that.
#
# Some examples of using the script:
# ./build_tests.sh -cicc -r -mrelease
# [The above line uses the intel compiler (icc) in a release config, and
# runs the tests after they're built.  No inline assembly is used in the build.]
#
# ./build_tests.sh -cclang -r -a
# [The above line uses the clang compiler in the default config (Debug), and
# runs the tests after they're built.  Inline assembly is used for the build
# wherever possible.]
#
# ./build_tests.sh -cgcc -mdebug
# [The above line uses the gcc compiler in a debug config.  Although it builds
# the tests, it does not run any of them.  No inline assembly is used.]



# ------------ How to get/install the various compilers on Linux --------------

# how to get gcc/g++ version 10, on ubuntu linux
# ----------------

# sudo add-apt-repository ppa:ubuntu-toolchain-r/test
# sudo apt-get update
# sudo apt install gcc-10
# sudo apt install g++-10
# sudo apt update
# sudo apt upgrade
# sudo apt install build-essential
# sudo apt update
# sudo apt upgrade
# sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 10
# sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 10
#  [ you may need to repeat the above two lines for any older gcc/g++ on your system ]
# sudo update-alternatives --config gcc
#  [choose the version of gcc you want as default]
# sudo update-alternatives --config g++
#  [choose the version of g++ you want as default]

#info taken from:
# https://askubuntu.com/questions/1192955/how-to-install-g-10-on-ubuntu-18-04
# https://www.fosslinux.com/39386/how-to-install-multiple-versions-of-gcc-and-g-on-ubuntu-20-04.htm
# https://linuxconfig.org/how-to-switch-between-multiple-gcc-and-g-compiler-versions-on-ubuntu-20-04-lts-focal-fossa



# how to get the latest version of clang (10 at the time of writing this)
# ----------------
# sudo apt update
# sudo apt upgrade
# bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)"
#   [The above line is what is recommended at https://apt.llvm.org/]
# sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-10 10
# sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-10 10
#  [ you may need to repeat the above two lines for any older clang/clang++ on
#    your system ]



# You get the Intel C++ compiler by downloading/installing Parallel Studio XE.
# It is free to install and use for open source projects.  See-
#   https://software.intel.com/content/www/us/en/develop/articles/qualify-for-free-software.html#opensourcecontributor
# After installing intel compiler on Ubuntu (I used Parallel Studio XE 2020
# Update1 Professional Edition), I did:
# ----------------
# sudo apt update
# sudo apt upgrade
# source /opt/intel/bin/compilervars.sh intel64
#   [ Whenever you wish to compile with Intel C++ (icc/icpc), you will need to
#     first run this compilervars script.  Probably you could set up your system
#     so that compilervars.sh runs at boot time and so that its effects apply to
#     any shell you start, but I would recommend *against* doing so.  The
#     compilervars script alters the PATH environment variable and puts some
#     intel bin directories as the first entries in the PATH.  Unfortunately
#     those bin directories contain clang, clang++, gcc, and g++ symbolic links
#     to (evidently) binaries that the Parallel Studio XE installer put on the
#     system.  This seems to me to be a bad choice by intel, because after
#     sourcing the compilervars script, the commands gcc, g++, clang, clang++ no
#     longer refer to any of the gcc or clang compilers that you may have
#     installed on your system.  Instead these commands will refer to the intel-
#     provided gcc/clang compilers, which will likely be a completely different
#     version of those compilers from what you expect or want.
#     The solution to this problem I am using is to start a subshell, then
#     source compilervars.sh, then perform the compilation with icc, then exit
#     back to the parent shell.  The commands I use to do this are:
#          bash
#          source /opt/intel/bin/compilervars.sh intel64
#          ... run icc, or run cmake with icc, etc ...
#          exit
#     Any running script is already in a subshell (unless you use source or . ),
#     so creating another subshell isn't necessary.  This present script here
#     assumes it is already in a subshell, and just sources compilervars.sh if
#     icc is chosen. ]
# sudo update-alternatives --config gcc
#   [ The above problems with compilervars.sh weren't my only issues.  For
#     whatever reason, I have experienced that the intel version of gcc/g++ is
#     not always the one used after I have sourced compilervars.sh, presumably
#     because the script altered the PATH variable to place the intel compiler
#     directories *after* usr/bin where the normal gcc/g++ on my system lies.
#     This is a problem for icc, because Intel C++ uses and needs gcc's system
#     includes and standard library: see  https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compatibility-and-portability/gcc-compatibility-and-interoperability.html
#     If icc mistakenly ends up using a version of gcc (for example, the very
#     latest version of gcc) that is not compatible with icc, then icc may fail
#     with compile errors that make it pretty obvious it is using your system
#     gcc rather than the intel-placed gcc.  Therefore, the default version of
#     gcc that is invoked by the commands 'gcc' and 'g++' *MUST NOT* be more
#     recent than your installation of Intel C++ supports.  I don't know what
#     the most recent gcc version is that icc from Parallel Studio XE 2020 (icc
#     v19.1) supports, but I know that when update-alternatives was set to
#     gcc-10 I had errors when I tried to compile icc with -std=c++17, and I
#     know that when I set update-alternatives to gcc-7 all the errors with
#     -std=c++17 went away.  Everything has worked fine for me with IntelC++
#     after I set the default 'gcc' and 'g++' to be version 7.  I would guess
#     that none of this is an issue if/when compilervars.sh places the intel
#     directories as the first entries in PATH, but I don't know.  I'm not sure
#     why compilervars.sh seems to be inconsistent about entry placement in PATH
#     for me.  Regardless, for whatever version of gcc that you are currently
#     using as the system default (via update-alteratives), you'll have to try
#     compiling with icc using -std=c++17 in order to know if your default gcc
#     version is making icc vulnerable to compile errors.  For example, with
#     this present script here when using icc, if you start getting unexpected
#     compile errors when you change this script's line that sets cpp_standard
#     to use cpp_standard="-std=c++17", then you'll know you need to default to
#     an older gcc version (via update-alternatives).  Either that, or figure
#     out how to make sure compilervars.sh always consistently places the intel
#     directories first in PATH. ]
# sudo update-alternatives --config g++
#   [ You *MUST* set the default g++ to be the same version as you just chose
#     for update-alternatives --config gcc.  For me, that meant g++-7. ]
#
# [ I didn't need to run "update-alternatives --config icc" (or icpc), but if
#   you install multiple versions of icc, then I expect you'd need to run
#   update-alternatives for both icc and icpc. ]



while getopts ":m:c:h-:raus" opt; do
  case $opt in
    h)
      ;&
    -)
      echo "Usage: build_tests [-c<compiler_name>] [-r] [-a] [-u] [-s] [-m<Release|Debug>]" >&2
      exit 1
      ;;
    c)
      compiler=$OPTARG
      ;;
    m)
      mode=$OPTARG
      ;;
    r)
      run_tests=true
      ;;
    a)
      use_inline_asm="-DHURCHALLA_ALLOW_INLINE_ASM_REDC=1"
      ;;
    u)
      use_all_inline_asm="-DHURCHALLA_ALLOW_INLINE_ASM_ALL=1"
      ;;
    s)
      test_avoid_cselect="-DHURCHALLA_AVOID_CSELECT=1"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z "$compiler" ]; then
  compiler=clang
fi
if [ -z "$mode" ]; then
  mode=Debug
fi


# We do heavy testing via this script.  Other methods for invoking the build and
# testing (e.g. github actions) usually aren't meant to be heavyweight and won't
# define this macro.
test_heavyweight="-DHURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT=1"


# Compiler commands
if [ "${compiler,,}" = "gcc" ] || [ "${compiler,,}" = "g++" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=g++
  cmake_c_compiler=-DCMAKE_C_COMPILER=gcc
  compiler_name=gcc
elif [ "${compiler,,}" = "gcc-7" ] || [ "${compiler,,}" = "g++-7" ] ||
     [ "${compiler,,}" = "gcc7" ] || [ "${compiler,,}" = "g++7" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=g++-7
  cmake_c_compiler=-DCMAKE_C_COMPILER=gcc-7
  compiler_name=gcc7
elif [ "${compiler,,}" = "gcc-10" ] || [ "${compiler,,}" = "g++-10" ] ||
     [ "${compiler,,}" = "gcc10" ] || [ "${compiler,,}" = "g++10" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=g++-10
  cmake_c_compiler=-DCMAKE_C_COMPILER=gcc-10
  compiler_name=gcc10
elif [ "${compiler,,}" = "clang" ] || [ "${compiler,,}" = "clang++" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=clang++
  cmake_c_compiler=-DCMAKE_C_COMPILER=clang
  compiler_name=clang
elif [ "${compiler,,}" = "clang-3" ] || [ "${compiler,,}" = "clang++-3" ] ||
     [ "${compiler,,}" = "clang3" ] || [ "${compiler,,}" = "clang++3" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=clang++-3.9
  cmake_c_compiler=-DCMAKE_C_COMPILER=clang-3.9
  compiler_name=clang3
elif [ "${compiler,,}" = "clang-6" ] || [ "${compiler,,}" = "clang++-6" ] ||
     [ "${compiler,,}" = "clang6" ] || [ "${compiler,,}" = "clang++6" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=clang++-6.0
  cmake_c_compiler=-DCMAKE_C_COMPILER=clang-6.0
  compiler_name=clang6
elif [ "${compiler,,}" = "clang-10" ] || [ "${compiler,,}" = "clang++-10" ] ||
     [ "${compiler,,}" = "clang10" ] || [ "${compiler,,}" = "clang++10" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=clang++-10
  cmake_c_compiler=-DCMAKE_C_COMPILER=clang-10
  compiler_name=clang10
elif [ "${compiler,,}" = "icc" ] || [ "${compiler,,}" = "icpc" ]; then
  cmake_cpp_compiler=-DCMAKE_CXX_COMPILER=icpc
  cmake_c_compiler=-DCMAKE_C_COMPILER=icc
  compiler_name=icc
  source /opt/intel/bin/compilervars.sh intel64
elif [ -n "$compiler" ]; then
  echo "Invalid argument for option -c: $compiler"
  exit 1
fi

echo Using compiler $compiler_name ...
echo Using build mode $mode ...



cpp_standard="-std=c++11"

# A long note about issues setting the C++ standard when using CMake
# ------------------------------------------------------------------
# In this bash file we're setting up the C++ standard as a compiler argument.
# This is not a normal or recommended method.  The normal/recommended method
# would be to set the standard in CMakeLists.txt.  However, in practice doing it
# that way is difficult to get to reliably work right for testing purposes, and
# more generally, it introduces potential issues with mixed standards.  For
# reference--- getting CMake to handle the standard would require two lines in
# CMakeLists.txt:
#    target_compile_features(MyTarget INTERFACE cxx_std_14)
#    set_target_properties(MyTarget PROPERTIES CXX_EXTENSIONS OFF)
# We'd need the second line since gcc will use -std=gnu++14 instead of
# -std=c++14, unless we turn off CXX_EXTENSIONS.  But these two things are still
# not enough in some circumstances: if we use cxx_std_11 instead of cxx_std_14
# above, cmake won't use any flags for the compiler standard for gcc (at least
# this is what I saw for gcc version 7) and this means that gcc will use gnu++11
# by default -- even though we specified cxx_std_11 and set CXX_EXTENSIONS off.
# It was at this point I gave up on setting the C++ standard within CMake.
#
# Instead, strictly for testing, I'm setting the standard as a compiler argument
# that gets passed to CMake.  Keep in mind that setting the standard via the
# command line is *only* done for testing.  My library CMakeLists.txt specifies
# nothing for the standard, and it will thus compile using whatever standard has
# been set (or has been defaulted) by a CMakeLists.txt that is using/consuming
# the library.  For these tests, the consuming CMakeLists.txt is
# test/CMakeLists.txt, but this CMakeLists.txt happens to be a special case of
# testing, where I need precise control of the standard and thus for the reasons
# given above I do it via compiler arguments.
#
# In principle it'd be easy to argue that my library CMakeLists.txt ought to
# specify that it needs at least C++11, and that it should do this via
# target_compile_features, but I have reasons above and below that go against
# doing that, and in practice in 2020 it would be fairly unusual for someone to
# use a C++ compiler (or to set the compiler standard in a way) that doesn't
# support C++11.  In such a case the worst thing that would happen, during
# compilation, is that the library would compile with errors that would probably
# be obvious missing C++11 features.  It'd much easier to make the argument that
# I ought to specify the standard in CMakeLists.txt for a library that needs
# C++17 or C++ 20, simply because it's much less likely (in 2020) that a
# compiler would happen to be set to use a high enough standard.  I'd still
# consider the related problem below though, even in that case.
# 
# A related problem is that setting the standard in my library's CMakeLists.txt
# could potentially lead to multiple standards being used to compile different
# parts of someone's project (e.g. C++11 for one library target and C++17 for
# another target, which then get linked together).  Although it'd be unusual for
# a problem to result from this, it isn't a good idea, and if/when it causes a
# problem it won't be obvious what went wrong.  So I avoid setting the standard
# anywhere in the CMakeLists.txt for the libraries I produce. There will still
# remain some scenarios where the standards could get mixed, but this reduces
# the exposure to it and it prevents my CMakeLists.txt from being a contributing
# cause of it.  For info on how/why mixing the standard could potentially result
# in errors, see:
# https://github.com/abseil/abseil-cpp/issues/259
# https://cgold.readthedocs.io/en/latest/tutorials/toolchain/globals/cxx-standard.html
# https://stackoverflow.com/questions/10717106/can-different-gcc-dialects-be-linked-together
# https://stackoverflow.com/questions/46746878/is-it-safe-to-link-c17-c14-and-c11-objects
# https://gcc.gnu.org/wiki/Cxx11AbiCompatibility
# https://cullmann.io/posts/cpp-standard-version-mix-up/




# Note: Gcc static analyzer was added in gcc version 10

#static analyzers
#-----
if [ "$compiler_name" = "gcc" ]; then
#  gcc_static_analysis="-fanalyzer"
# !!! g++-10 with -fanalyzer locks the system on gtest !!!
# g++-10 with -Wanalyzer-too-complex had compiler errors on gtest

  : # do nothing, at least for now

elif [ "$compiler_name" = "clang" ]; then
#  clang_static_analysis=(-DCMAKE_CXX_CLANG_TIDY="clang-tidy;-checks=-*,clang-analyzer-*")
  : # do nothing, at least for now
fi


#undefined behavior sanitizers
#-----
if [ "$compiler_name" = "gcc" ]; then
  gcc_ubsan="-fsanitize=undefined -fno-sanitize-recover \
           -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow"

elif [ "$compiler_name" = "clang" ]; then
  # My installed version of clang doesn't support -fsanitize=implicit-conversion
  clang_ubsan="-fsanitize=undefined -fsanitize=nullability -fsanitize=bounds \
             -fsanitize=float-divide-by-zero"

  # The next line in a perfect world wouldn't be needed, but for some versions
  # of clang (clang 10 for me), the linker doesn't find __muloti4 when using the
  # undefined behavior sanitizers.  __muloti4 is defined in compiler-rt.
  # See https://bugs.llvm.org/show_bug.cgi?id=16404
  clang_ubsan_link_flags="-rtlib=compiler-rt -lgcc_s"
fi


#address sanitizers
#-----
clang_asan=""
gcc_asan="-fsanitize=address"
# clang -fsanitize=address -O1 -fno-omit-frame-pointer -g   tests/use-after-free.c
# clang++ -O1 -g -fsanitize=address -fno-omit-frame-pointer example_UseAfterFree.cc

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# but modern cmake would be to use target_compile_options

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined")

#target_link_libraries(MyTarget
#  -fsanitize=address
#)

# prefer target_link_options instead?

#target_compile_options(your_target
#  PRIVATE
#  $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Debug>:${CMAKE_CXX_FLAGS_RELEASE}>>

# target_compile_options(tgt PRIVATE "/MD$<$<CONFIG:Debug>:d>")


# target_compile_options(foo PUBLIC "$<$<CONFIG:DEBUG>:${MY_DEBUG_OPTIONS}>")
# target_compile_options(foo PUBLIC "$<$<CONFIG:RELEASE>:${MY_RELEASE_OPTIONS}>")

# target_compile_options(MyLib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-std=c++14>)

# target_compile_options(MyLib PUBLIC "$<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:DEBUG>>:/MDd>")


# if(MSVC)
#     add_compile_options("/W4" "$<$<CONFIG:RELEASE>:/O2>")
# else()
#     add_compile_options("-Wall" "-Wextra" "-Werror" "$<$<CONFIG:RELEASE>:-O3>")
#     if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#         add_compile_options("-stdlib=libc++")
#     else()
#         # nothing special for gcc at the moment
#     endif()
# endif()



#modes
# 1. Asan+UBsan+Lsan
# 2. Tsan
# 3. Msan
# 4. Valgrind (can't be used with other sanitizers)

# a run of "splint" and/or cppcheck
# cpplint
# include what you use (iwyu), and lwyu
# Clang-Tidy
# CppCoreCheck

# doxygen

# <LANG>_CLANG_TIDY: CMake 3.6+
# <LANG>_CPPCHECK
# <LANG>_CPPLINT
# <LANG>_INCLUDE_WHAT_YOU_USE
# LINK_WHAT_YOU_USE




exit_on_failure () {
  if [ $? -ne 0 ]; then
    exit 1
  fi
}


# https://stackoverflow.com/questions/59895/how-to-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ "${mode,,}" = "release" ]; then
    pushd script_dir > /dev/null 2>&1
    build_dir=build/release_$compiler_name
    mkdir -p $build_dir
    cmake -S. -B./$build_dir -DTEST_HURCHALLA_LIBS=ON \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_CXX_FLAGS="$cpp_standard  \
            $test_avoid_cselect  $test_heavyweight \
            $use_inline_asm  $use_all_inline_asm \
            $gcc_static_analysis"  "${clang_static_analysis[@]}" \
            $cmake_cpp_compiler $cmake_c_compiler
    exit_on_failure
    cmake --build ./$build_dir --config Release
    exit_on_failure
    popd > /dev/null 2>&1
elif [ "${mode,,}" = "debug" ]; then
    pushd script_dir > /dev/null 2>&1
    build_dir=build/debug_$compiler_name
    mkdir -p $build_dir
    cmake -S. -B./$build_dir -DTEST_HURCHALLA_LIBS=ON \
            -DCMAKE_BUILD_TYPE=Debug \
            -DCMAKE_EXE_LINKER_FLAGS="$clang_ubsan_link_flags" \
            -DCMAKE_CXX_FLAGS="$cpp_standard  $clang_ubsan  $gcc_ubsan  \
            $test_avoid_cselect  $test_heavyweight \
            $use_inline_asm  $use_all_inline_asm \
            $gcc_static_analysis"  "${clang_static_analysis[@]}" \
            $cmake_cpp_compiler $cmake_c_compiler
    exit_on_failure
    cmake --build ./$build_dir --config Debug
    exit_on_failure
    popd > /dev/null 2>&1
else
    echo "Invalid argument for option -m: $mode"
    exit 1
fi


# -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
# cmake  -S.  -B./build_tmp  -DCMAKE_CXX_FLAGS="-std=c++17"  -DTEST_HURCHALLA_LIBS=ON  -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_CXX_COMPILER=icpc  -DCMAKE_C_COMPILER=icc
# cmake --build ./build_tmp --config Debug


if [ "$run_tests" = true ]; then
  ./$build_dir/test_ndebug_programming_by_contract --gtest_break_on_failure
  exit_on_failure
  ./$build_dir/test_programming_by_contract --gtest_break_on_failure
  exit_on_failure
  ./$build_dir/test_hurchalla_util --gtest_break_on_failure
  exit_on_failure
  ./$build_dir/test_hurchalla_modular_arithmetic --gtest_break_on_failure
  exit_on_failure
fi
