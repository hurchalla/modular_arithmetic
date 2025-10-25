#!/bin/bash

# Copyright (c) 2025 Jeffrey Hurchalla.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


invoke_test() {
  (
    set -e
    echo ""
    shift
    echo $@
    echo ""
    $@
  ) 2>&1 | tee $1

  # Check if the subshell failed
  if [ "${PIPESTATUS[0]}" -ne 0 ]; then
    echo "One of the commands failed. Exiting."
    exit 1
  fi
}



invoke_test 64_quarter_gcc_asm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_quarter_gcc_asm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_quarter_gcc_asm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_quarter_clang_asm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_quarter_clang_asm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_quarter_clang_asm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_quarter_gcc_noasm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_SCALAR

invoke_test 64_quarter_gcc_noasm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_ARRAY

invoke_test 64_quarter_gcc_noasm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY

invoke_test 64_quarter_clang_noasm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_SCALAR

invoke_test 64_quarter_clang_noasm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_ARRAY

invoke_test 64_quarter_clang_noasm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY



invoke_test 64_half_gcc_asm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_half_gcc_asm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_half_gcc_asm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_half_clang_asm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_half_clang_asm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_half_clang_asm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_half_gcc_noasm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_SCALAR

invoke_test 64_half_gcc_noasm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_ARRAY

invoke_test 64_half_gcc_noasm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY

invoke_test 64_half_clang_noasm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_SCALAR

invoke_test 64_half_clang_noasm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_ARRAY

invoke_test 64_half_clang_noasm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY



invoke_test 64_full_gcc_asm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_full_gcc_asm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_full_gcc_asm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_full_clang_asm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_full_clang_asm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_full_clang_asm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 64_full_gcc_noasm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_SCALAR

invoke_test 64_full_gcc_noasm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_ARRAY

invoke_test 64_full_gcc_noasm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY

invoke_test 64_full_clang_noasm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_SCALAR

invoke_test 64_full_clang_noasm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_ARRAY

invoke_test 64_full_clang_noasm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull uint64_t 191 8 22 \
      -DTEST_PARTIAL_ARRAY






invoke_test 128_quarter_gcc_asm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_quarter_gcc_asm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_quarter_gcc_asm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_quarter_clang_asm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_quarter_clang_asm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_quarter_clang_asm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_quarter_gcc_noasm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_SCALAR

invoke_test 128_quarter_gcc_noasm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_ARRAY

invoke_test 128_quarter_gcc_noasm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY

invoke_test 128_quarter_clang_noasm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_SCALAR

invoke_test 128_quarter_clang_noasm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_ARRAY

invoke_test 128_quarter_clang_noasm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryQuarter __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY



invoke_test 128_half_gcc_asm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_half_gcc_asm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_half_gcc_asm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_half_clang_asm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_half_clang_asm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_half_clang_asm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_half_gcc_noasm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_SCALAR

invoke_test 128_half_gcc_noasm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_ARRAY

invoke_test 128_half_gcc_noasm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY

invoke_test 128_half_clang_noasm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_SCALAR

invoke_test 128_half_clang_noasm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_ARRAY

invoke_test 128_half_clang_noasm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryHalf __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY



invoke_test 128_full_gcc_asm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_full_gcc_asm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_full_gcc_asm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_full_clang_asm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_SCALAR -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_full_clang_asm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_full_clang_asm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY -DHURCHALLA_MONTGOMERY_TWO_POW_USE_CSELECT_ON_BIT -DHURCHALLA_ALLOW_INLINE_ASM_ALL

invoke_test 128_full_gcc_noasm_scalar.txt ./testbench_2kary.sh g++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_SCALAR

invoke_test 128_full_gcc_noasm_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_ARRAY

invoke_test 128_full_gcc_noasm_partial_array.txt ./testbench_2kary.sh g++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY

invoke_test 128_full_clang_noasm_scalar.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_SCALAR

invoke_test 128_full_clang_noasm_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_ARRAY

invoke_test 128_full_clang_noasm_partial_array.txt ./testbench_2kary.sh clang++ O3 MontgomeryFull __uint128_t 191 8 40 \
      -DTEST_PARTIAL_ARRAY
