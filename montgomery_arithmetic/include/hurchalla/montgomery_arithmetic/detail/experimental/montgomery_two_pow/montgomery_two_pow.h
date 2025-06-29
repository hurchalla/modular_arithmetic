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
// Implementation note: uses a modification of the k-ary exponentiation alg
// ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring )
// which precalculates the even exponents as well as the normal odd exponents,
// in order to avoid two conditional branches that would exist in the main loop
// of the normal k-ary algorithm.
//
template <class MF, typename U>
typename MF::MontgomeryValue montgomery_two_pow(const MF& mf, U n)
{
    static_assert(ut_numeric_limits<U>::is_integer);
    static_assert(!ut_numeric_limits<U>::is_signed);

    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using std::size_t;

    constexpr int P = 5;   // (1 << P) == the k in k-ary exponentiation

    // initialize the precalculation table for k-ary pow algorithm
    static_assert(P > 0);
    constexpr size_t TABLESIZE = 1 << P;   
    C table[TABLESIZE];
    table[0] = mf.getUnityValue();   // montgomery one
    for (size_t i=1; i<TABLESIZE; ++i)
        table[i] = mf.add(table[i-1], table[i-1]);

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

    int shift = numbits - P;
    U tmp = n >> shift;
    HPBC_ASSERT(tmp <= MASK);
    // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
    size_t index = static_cast<size_t>(tmp);
    V result = table[index];


    while (shift >= P) {
//        HURCHALLA_REQUEST_UNROLL_LOOP
        for (int i=0; i<P; ++i)
            result = mf.square(result);
#if 1
// sliding window optimization
// TODO: maybe optimize the shift in the next line, since n may be > 64bit
        while (shift > P && (static_cast<size_t>(n >> (shift-1)) & 1) == 0) {
            result = mf.square(result);
            --shift;
        }
#endif
        shift -= P;
// TODO: maybe optimize next line, since n may be > 64bit
        tmp = n >> shift;
        index = static_cast<size_t>(tmp) & MASK;
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




template <class MF, typename U, size_t ARRAY_SIZE>
std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
array_montgomery_two_pow(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
{
    static_assert(ut_numeric_limits<U>::is_integer);
    static_assert(!ut_numeric_limits<U>::is_signed);

    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using std::size_t;

    constexpr int P = 5;   // (1 << P) == the k in k-ary exponentiation

    // initialize the precalculation table for k-ary pow algorithm
    static_assert(P > 0);
    constexpr size_t TABLESIZE = 1 << P;   
    C table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
        table[0][j] = mf[j].getUnityValue();   // montgomery one
    for (size_t i=1; i<TABLESIZE; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            table[i][j] = mf[j].add(table[i-1][j], table[i-1][j]);
    }

    U n_max = n[0];
    HURCHALLA_REQUEST_UNROLL_LOOP for (int j=1; j<ARRAY_SIZE; ++j) {
        if (n_max < n[j])
            n_max = n[j];
    }

    constexpr size_t MASK = TABLESIZE - 1;
    if (n_max <= MASK) {
        std::array<V, ARRAY_SIZE> result;
        HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
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

    int shift = numbits - P;
    std::array<V, ARRAY_SIZE> result;
    std::array<U, ARRAY_SIZE> tmp;
    HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
        tmp[j] = n[j] >> shift;
        HPBC_ASSERT(tmp[j] <= MASK);
        // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
        result[j] = table[static_cast<size_t>(tmp[j])][j];
    }


    while (shift >= P) {
        for (int i=0; i<P; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].square(result[j]);
        }
        shift -= P;
        std::array<size_t, ARRAY_SIZE> index;
        HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
// TODO: maybe optimize next line, since n may be > 64bit
            tmp[j] = n[j] >> shift;
            index[j] = static_cast<size_t>(tmp[j]) & MASK;
            result[j] = mf[j].multiply(result[j], table[index[j]][j]);
        }
    }


    if (shift == 0)
        return result;
    HPBC_ASSERT(0 < shift && shift < P);

    for (int i=0; i<shift; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            result[j] = mf[j].square(result[j]);
    }
    size_t tmpmask = (1 << shift) - 1;
    std::array<size_t, ARRAY_SIZE> index;
    HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
        index[j] = static_cast<size_t>(n[j]) & tmpmask;
        result[j] = mf[j].multiply(result[j], table[index[j]][j]);
    }
    return result;
}




} // end namespace

#endif
