// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_TWO_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_TWO_POW_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/count_leading_zeros.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>
#include <cstddef>
#include <array>

namespace hurchalla {


// Calculates pow(2, n), modulo the modulus of mf, and returns the result in
// montgomeryform representation.
//
// Implementation note: this is a modified version of the 2^k-ary exponentiation
// algorithm  ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring ),
// which precalculates the even exponents as well as the normal odd exponents,
// in order to avoid two conditional branches that would exist in the main loop
// of the normal 2^k-ary algorithm.  This is particularly helpful for the array
// version of this function further below.
//
template <class MF, typename U,
          bool USE_SLIDING_WINDOW_OPTIMIZATION = true,
          size_t TABLE_BITS = 5>
typename MF::MontgomeryValue montgomery_two_pow(const MF& mf, U n)
{
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    // FYI you almost certainly want 2 <= TABLE_BITS <= 6.
    // TABLE_BITS > 0 is required; anything above 9 is probably a very bad idea
    // even if it works (it would mean calculating 1024+ table entries!)
    static_assert(0 < TABLE_BITS && TABLE_BITS < 10, "");

    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using std::size_t;

    constexpr size_t P = TABLE_BITS;

    // initialize the precalculation table for 2^k-ary pow algorithm
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1 << P;
    V table[TABLESIZE];
    static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
    table[0] = mf.getUnityValue();   // montgomery one

    // This simple 'for' loop would work fine, but we can init faster...
//    for (int i=1; i<TABLESIZE; ++i)
//        table[i] = mf.add(table[i-1], table[i-1]);

    // Let's do optimized initializations for the different table sizes.
    if (TABLESIZE == 2) {
#if 0
        table[1] = mf.add(table[0], table[0]);
#else
// For this particular case of TABLE_BITS == 1 (TABLESIZE == 2), we can
// use a version of 2^k-ary that is heavily optimized for the 1 bit table:
        if (n <= 1) {
            if (n == 0)
                return mf.getUnityValue();
            else
                return mf.add(mf.getUnityValue(), mf.getUnityValue());
        }
        HPBC_ASSERT(n > 1);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_ASSERT(numbits > 1);
        int shift = numbits - 1;
        HPBC_ASSERT((n >> shift) == 1);   // because shift == numbits - 1
        C mont_two = mf.add(mf.getUnityValue(), mf.getUnityValue());
        V result = mont_two;
        HPBC_ASSERT(shift >= 1);
        result = mf.add(result, result);  // since result == 2, this add is equivalent to squaring
        --shift;
        V vtmp = mf.add(result, result);
        result.cmov(static_cast<size_t>(n >> shift) & 1, vtmp);

        while (shift >= 1) {
            result = mf.square(result);
            --shift;
            vtmp = mf.add(result, result);
            result.cmov(static_cast<size_t>(n >> shift) & 1, vtmp);
        }
        return result;
#endif
    }
    if (TABLESIZE == 4) {
        table[1] = mf.add(table[0], table[0]);
        table[2] = mf.add(table[1], table[1]);
        table[3] = mf.add(table[2], table[2]);
    } else if (TABLESIZE == 8) {
        table[1] = mf.add(table[0], table[0]);
        table[2] = mf.add(table[1], table[1]);
        table[3] = mf.add(table[2], table[2]);
        table[4] = mf.add(table[3], table[3]);
        table[5] = mf.add(table[4], table[4]);
        table[6] = mf.add(table[5], table[5]);
        table[7] = mf.add(table[6], table[6]);
    } else if (TABLESIZE == 16) {
        table[1] = mf.add(table[0], table[0]);
        table[2] = mf.add(table[1], table[1]);
        table[3] = mf.add(table[2], table[2]);
        table[4] = mf.add(table[3], table[3]);
        table[5] = mf.add(table[4], table[4]);
        table[6] = mf.add(table[5], table[5]);
        table[7] = mf.add(table[6], table[6]);

        table[14] = mf.square(table[7]);

        table[8]  = mf.add(table[7], table[7]);
        table[9]  = mf.add(table[8], table[8]);
        table[10] = mf.add(table[9], table[9]);
        table[11] = mf.add(table[10], table[10]);
        table[12] = mf.add(table[11], table[11]);
        table[13] = mf.add(table[12], table[12]);

        table[15] = mf.add(table[14], table[14]);
    }
    else if (TABLESIZE == 32) {
// these #if #elifs should be functionally equivalent.  You can test to see
// which is fastest.
#if 0
        for (size_t i=0; i<31; ++i)
            table[i+1] = mf.add(table[i], table[i]);
#elif 0
        // Calling convertIn() is dubious- it ties us to using a MontgomeryForm
        // that computes RSquaredModN in its constructor, which is intensive and
        // something that in principle we can avoid if we never call convertIn.
        table[14] = mf.convertIn(64 * 256);
        table[23] = mf.convertIn(128 * 256 * 256);

        table[1] = mf.add(table[0], table[0]);
        table[2] = mf.add(table[1], table[1]);
        table[3] = mf.add(table[2], table[2]);
        table[4] = mf.add(table[3], table[3]);
        table[5] = mf.add(table[4], table[4]);

        table[6] = mf.add(table[5], table[5]);
        table[15] = mf.add(table[14], table[14]);
        table[24] = mf.add(table[23], table[23]);

        table[7] = mf.add(table[6], table[6]);
        table[16] = mf.add(table[15], table[15]);
        table[25] = mf.add(table[24], table[24]);

        table[8]  = mf.add(table[7], table[7]);
        table[17] = mf.add(table[16], table[16]);
        table[26] = mf.add(table[25], table[25]);

        table[9]  = mf.add(table[8], table[8]);
        table[18] = mf.add(table[17], table[17]);
        table[27] = mf.add(table[26], table[26]);

        table[10] = mf.add(table[9], table[9]);
        table[19] = mf.add(table[18], table[18]);
        table[28] = mf.add(table[27], table[27]);

        table[11] = mf.add(table[10], table[10]);
        table[20] = mf.add(table[19], table[19]);
        table[29] = mf.add(table[28], table[28]);

        table[12] = mf.add(table[11], table[11]);
        table[21] = mf.add(table[20], table[20]);
        table[30] = mf.add(table[29], table[29]);

        table[13] = mf.add(table[12], table[12]);
        table[22] = mf.add(table[21], table[21]);
        table[31] = mf.add(table[30], table[30]);
#elif 1
// This tested fastest for me when using __uint128_t
        table[1] = mf.add(table[0], table[0]);
        table[2] = mf.add(table[1], table[1]);
        table[3] = mf.add(table[2], table[2]);
        table[4] = mf.add(table[3], table[3]);
        table[5] = mf.add(table[4], table[4]);
        table[6] = mf.add(table[5], table[5]);
        table[7] = mf.add(table[6], table[6]);
        table[8]  = mf.add(table[7], table[7]);
        table[9]  = mf.add(table[8], table[8]);
        table[10] = mf.add(table[9], table[9]);
        table[11] = mf.add(table[10], table[10]);
        table[12] = mf.add(table[11], table[11]);
        table[13] = mf.add(table[12], table[12]);

        table[26] = mf.square(table[13]);
//        table[26] = mf.convertIn(4 * 256 * 256 * 256);
        table[14] = mf.add(table[13], table[13]);
        table[15] = mf.add(table[14], table[14]);
        table[16] = mf.add(table[15], table[15]);
        table[17] = mf.add(table[16], table[16]);
        table[18] = mf.add(table[17], table[17]);
        table[19] = mf.add(table[18], table[18]);
        table[20] = mf.add(table[19], table[19]);

        table[21] = mf.add(table[20], table[20]);
        table[27] = mf.add(table[26], table[26]);

        table[22] = mf.add(table[21], table[21]);
        table[28] = mf.add(table[27], table[27]);

        table[23] = mf.add(table[22], table[22]);
        table[29] = mf.add(table[28], table[28]);

        table[24] = mf.add(table[23], table[23]);
        table[30] = mf.add(table[29], table[29]);

        table[25] = mf.add(table[24], table[24]);
        table[31] = mf.add(table[30], table[30]);
#elif 1
        table[1] = mf.add(table[0], table[0]);
        table[2] = mf.add(table[1], table[1]);
        table[3] = mf.add(table[2], table[2]);
        table[4] = mf.add(table[3], table[3]);
        table[5] = mf.add(table[4], table[4]);
        table[6] = mf.add(table[5], table[5]);
        table[7] = mf.add(table[6], table[6]);
        table[8]  = mf.add(table[7], table[7]);
        table[9]  = mf.add(table[8], table[8]);

        table[18] = mf.square(table[9]);
        table[10] = mf.add(table[9], table[9]);
        table[11] = mf.add(table[10], table[10]);
        table[12] = mf.add(table[11], table[11]);
        table[13] = mf.add(table[12], table[12]);

        table[27] = mf.multiply(table[18], table[9]);
        table[14] = mf.add(table[13], table[13]);
        table[19] = mf.add(table[18], table[18]);
        table[15] = mf.add(table[14], table[14]);
        table[20] = mf.add(table[19], table[19]);
        table[16] = mf.add(table[15], table[15]);
        table[21] = mf.add(table[20], table[20]);
        table[17] = mf.add(table[16], table[16]);
        table[22] = mf.add(table[21], table[21]);

        table[23] = mf.add(table[22], table[22]);
        table[28] = mf.add(table[27], table[27]);
        table[24] = mf.add(table[23], table[23]);
        table[29] = mf.add(table[28], table[28]);
        table[25] = mf.add(table[24], table[24]);
        table[30] = mf.add(table[29], table[29]);
        table[26] = mf.add(table[25], table[25]);
        table[31] = mf.add(table[30], table[30]);
#else
        table[1] = mf.add(table[0], table[0]);
        table[2] = mf.add(table[1], table[1]);
        table[3] = mf.add(table[2], table[2]);
        table[4] = mf.add(table[3], table[3]);
        table[5] = mf.add(table[4], table[4]);
        table[6] = mf.add(table[5], table[5]);
        table[7] = mf.add(table[6], table[6]);
        table[8] = mf.add(table[7], table[7]);

        table[16] = mf.square(table[8]);
        table[9]  = mf.add(table[8], table[8]);
        table[10] = mf.add(table[9], table[9]);
        table[11] = mf.add(table[10], table[10]);
        table[12] = mf.add(table[11], table[11]);

        table[24] = mf.square(table[12]);
        table[28] = mf.multiply(table[12], table[16]);

        table[13] = mf.add(table[12], table[12]);
        table[17] = mf.add(table[16], table[16]);
        table[14] = mf.add(table[13], table[13]);
        table[18] = mf.add(table[17], table[17]);
        table[15] = mf.add(table[14], table[14]);
        table[19] = mf.add(table[18], table[18]);

        table[20] = mf.add(table[19], table[19]);

        table[21] = mf.add(table[20], table[20]);
        table[25] = mf.add(table[24], table[24]);
        table[29] = mf.add(table[28], table[28]);

        table[22] = mf.add(table[21], table[21]);
        table[26] = mf.add(table[25], table[25]);
        table[30] = mf.add(table[29], table[29]);

        table[23] = mf.add(table[22], table[22]);
        table[27] = mf.add(table[26], table[26]);
        table[31] = mf.add(table[30], table[30]);
#endif
    }
    else if (TABLESIZE == 64) {
#if 0
        for (size_t i=0; i<63; ++i)
            table[i+1] = mf.add(table[i], table[i]);
#elif 0
        for (size_t i=0; i<23; ++i)
            table[i+1] = mf.add(table[i], table[i]);
        table[46] = mf.square(table[23]);
        table[24] = mf.add(table[23], table[23]);
        table[25] = mf.add(table[24], table[24]);
        table[26] = mf.add(table[25], table[25]);
        table[27] = mf.add(table[26], table[26]);
        table[28] = mf.add(table[27], table[27]);
        for (size_t i=28; i<45; ++i) {
            table[i+1] = mf.add(table[i], table[i]);
            table[i+19] = mf.add(table[i+18], table[i+18]);
        }
#else
// This tested fastest for me when using __uint128_t
        for (size_t i=0; i<16; ++i)
            table[i+1] = mf.add(table[i], table[i]);

        table[32] = mf.square(table[16]);

        table[17] = mf.add(table[16], table[16]);
        table[18] = mf.add(table[17], table[17]);
        table[19] = mf.add(table[18], table[18]);
        table[20] = mf.add(table[19], table[19]);
        table[21] = mf.add(table[20], table[20]);

        table[22] = mf.add(table[21], table[21]);
        table[33] = mf.add(table[32], table[32]);
        table[53] = mf.multiply(table[32], table[21]);

        for (size_t i=22; i<31; ++i) {
            table[i+1] = mf.add(table[i], table[i]);
            table[i+12] = mf.add(table[i+11], table[i+11]);
        }
        for (size_t i=42; i<52; ++i) {
            table[i+1] = mf.add(table[i], table[i]);
            table[i+12] = mf.add(table[i+11], table[i+11]);
        }
#endif
    }
    else {
        for (size_t i=0; i<TABLESIZE-1; ++i)
            table[i+1] = mf.add(table[i], table[i]);
    }


    constexpr size_t MASK = TABLESIZE - 1;
    if (n <= MASK)
        return table[static_cast<size_t>(n)];


    // count_leading_zeros returns the number of leading 0-bits in n, starting
    // at the most significant bit position. If n is 0, the result is undefined
    HPBC_ASSERT(n > 0);
    int leading_zeros = count_leading_zeros(n);
    int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
    // because we returned above if (n <= MASK), we can assert the following:
    HPBC_ASSERT(numbits > P);

    int shift = numbits - static_cast<int>(P);
    U tmp = n >> shift;
    HPBC_ASSERT(tmp <= MASK);
    // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
    size_t index = static_cast<size_t>(tmp);
    V result = table[index];


    while (shift >= static_cast<int>(P)) {
        if (USE_SLIDING_WINDOW_OPTIMIZATION) {
            while (shift > P && (static_cast<size_t>(n >> (shift-1)) & 1) == 0) {
                result = mf.square(result);
                --shift;
            }
        }

//        HURCHALLA_REQUEST_UNROLL_LOOP
        for (size_t i=0; i<P; ++i)
            result = mf.square(result);

        shift -= static_cast<int>(P);
        index = static_cast<size_t>(n >> shift) & MASK;
        result = mf.multiply(result, table[index]);
    }


    if (shift == 0)
        return result;
    HPBC_ASSERT(0 < shift && shift < P);

    for (int i=0; i<shift; ++i)
        result = mf.square(result);
    size_t tmpmask = (1 << shift) - 1;
    index = static_cast<size_t>(n) & tmpmask;
    result = mf.multiply(result, table[index]);
    return result;
}




template <class MF, typename U,
          size_t ARRAY_SIZE,
          size_t TABLE_BITS = 5>
std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
array_montgomery_two_pow(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
{
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    // FYI you almost certainly want 2 <= TABLE_BITS <= 6.
    // TABLE_BITS > 0 is required; anything above 9 is probably a very bad idea
    // even if it works (it would mean calculating 1024+ table entries!)
    static_assert(0 < TABLE_BITS && TABLE_BITS < 10, "");

    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using std::size_t;

    constexpr size_t P = TABLE_BITS;

    // Initialize the precalculation table...
    // We'll likely get very good instruction level parallelism via the array,
    // and have no need for any of the init tricks seen in the function above.
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1 << P;   
    C table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
        table[0][j] = mf[j].getUnityValue();   // montgomery one
    for (int i=1; i<TABLESIZE; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            table[i][j] = mf[j].add(table[i-1][j], table[i-1][j]);
    }

    U n_max = n[0];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=1; j<ARRAY_SIZE; ++j) {
        if (n_max < n[j])
            n_max = n[j];
    }

    constexpr size_t MASK = TABLESIZE - 1;
    if (n_max <= MASK) {
        std::array<V, ARRAY_SIZE> result;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = table[static_cast<size_t>(n[j])][j];
        return result;
    }


    // count_leading_zeros returns the number of leading 0-bits in n_max, starting
    // at the most significant bit position. If n_max is 0, the result is undefined
    HPBC_ASSERT(n_max > 0);
    int leading_zeros = count_leading_zeros(n_max);
    int numbits = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
    // because we returned above if (n_max <= MASK), we can assert the following:
    HPBC_ASSERT(numbits > P);

    int shift = numbits - static_cast<int>(P);
    std::array<V, ARRAY_SIZE> result;
    std::array<U, ARRAY_SIZE> tmp;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        tmp[j] = n[j] >> shift;
        HPBC_ASSERT(tmp[j] <= MASK);
        // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
        result[j] = table[static_cast<size_t>(tmp[j])][j];
    }


    while (shift >= static_cast<int>(P)) {
        for (size_t i=0; i<P; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<LowuopsTag>(result[j]);
        }
        shift -= static_cast<int>(P);
        std::array<size_t, ARRAY_SIZE> index;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            tmp[j] = n[j] >> shift;
            index[j] = static_cast<size_t>(tmp[j]) & MASK;
            result[j] = mf[j].template multiply<LowuopsTag>(result[j], table[index[j]][j]);
        }
    }


    if (shift == 0)
        return result;
    HPBC_ASSERT(0 < shift && shift < P);

    for (int i=0; i<shift; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mf[j].template square<LowuopsTag>(result[j]);
    }
    size_t tmpmask = (1 << shift) - 1;
    std::array<size_t, ARRAY_SIZE> index;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        index[j] = static_cast<size_t>(n[j]) & tmpmask;
        result[j] = mf[j].template multiply<LowuopsTag>(result[j], table[index[j]][j]);
    }
    return result;
}




} // end namespace

#endif
