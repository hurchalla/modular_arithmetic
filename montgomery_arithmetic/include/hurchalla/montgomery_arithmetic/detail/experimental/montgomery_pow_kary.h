// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_KARY_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_KARY_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/count_leading_zeros.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include <type_traits>
#include <cstddef>
#include <array>

namespace hurchalla {


// MF should be a MontgomeryForm type (i.e. either plain MontgomeryForm, or one
//   of its aliases like MontgomeryQuarter or MontgomeryHalf).  The implicit
//   modulus that this function uses is therefore equal to mf.getModulus().
// x is the base (in Montgomery domain)
// nexp is the exponent (it can be any integer type T).
//
// This function returns the Montgomery domain result of  x^nexp (mod modulus).
//
// Implementation note: uses a modification of the k-ary exponentiation alg
// ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring )
// which precalculates the even exponents as well as the normal odd exponents,
// in order to avoid two conditional branches that would exist in the main loop
// of the normal k-ary algorithm.
//
template <class MF, typename T>
typename MF::MontgomeryValue montgomery_pow_kary(const MF& mf, typename MF::MontgomeryValue x, T nexp)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    using U = typename extensible_make_unsigned<T>::type;

    HPBC_PRECONDITION(nexp >= 0);
    U n = static_cast<U>(nexp);

    using V = typename MF::MontgomeryValue;
    using std::size_t;

    constexpr size_t P = 4;   // (1 << P) == the k in k-ary exponentiation

    // initialize the precalculation table for k-ary pow algorithm
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1 << P;   
    static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
    V table[TABLESIZE];
    table[0] = mf.getUnityValue();   // montgomery one
    table[1] = x;
    if (TABLESIZE == 4) {
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
    } else if (TABLESIZE == 8) {
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
        table[4] = mf.square(table[2]);
        table[5] = mf.multiply(table[2], table[3]);
        table[6] = mf.square(table[3]);
        table[7] = mf.multiply(table[3], table[4]);
    } else if (TABLESIZE == 16) {
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
        table[4] = mf.square(table[2]);
        table[5] = mf.multiply(table[2], table[3]);
        table[6] = mf.square(table[3]);
        table[7] = mf.multiply(table[3], table[4]);
        table[8] = mf.square(table[4]);
        table[9] = mf.multiply(table[4], table[5]);
        table[10] = mf.square(table[5]);
        table[11] = mf.multiply(table[5], table[6]);
        table[12] = mf.square(table[6]);
        table[13] = mf.multiply(table[6], table[7]);
        table[14] = mf.square(table[7]);
        table[15] = mf.multiply(table[7], table[8]);
    } else if (TABLESIZE == 32) {
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
        table[4] = mf.square(table[2]);
        table[5] = mf.multiply(table[2], table[3]);
        table[6] = mf.square(table[3]);
        table[7] = mf.multiply(table[3], table[4]);
        table[8] = mf.square(table[4]);
        table[9] = mf.multiply(table[4], table[5]);
        table[10] = mf.square(table[5]);
        table[11] = mf.multiply(table[5], table[6]);
        table[12] = mf.square(table[6]);
        table[13] = mf.multiply(table[6], table[7]);
        table[14] = mf.square(table[7]);
        table[15] = mf.multiply(table[7], table[8]);

        table[16] = mf.square(table[8]);
        table[17] = mf.multiply(table[8], table[9]);
        table[18] = mf.square(table[9]);
        table[19] = mf.multiply(table[9], table[10]);
        table[20] = mf.square(table[10]);
        table[21] = mf.multiply(table[10], table[11]);
        table[22] = mf.square(table[11]);
        table[23] = mf.multiply(table[11], table[12]);
        table[24] = mf.square(table[12]);
        table[25] = mf.multiply(table[12], table[13]);
        table[26] = mf.square(table[13]);
        table[27] = mf.multiply(table[13], table[14]);
        table[28] = mf.square(table[14]);
        table[29] = mf.multiply(table[14], table[15]);
        table[30] = mf.square(table[15]);
        table[31] = mf.multiply(table[15], table[16]);
    } else {
        // we should check for a ((power of 2) >= 64), but this is
        // probably adquate or our needs
        static_assert(TABLESIZE % 64 == 0, "");
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=4; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            table[i] = mf.square(table[halfi]);
            table[i+1] = mf.multiply(table[halfi], table[halfi + 1]);
        }
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
//        HURCHALLA_REQUEST_UNROLL_LOOP
        for (size_t i=0; i<P; ++i)
            result = mf.square(result);
#if 1
// sliding window optimization
// TODO: maybe optimize the shift in the next line, since n may be > 64bit
        while (shift > P && (static_cast<size_t>(n >> (shift-1)) & 1) == 0) {
            result = mf.square(result);
            --shift;
        }
#endif
        shift -= static_cast<int>(P);
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






// This is the array version of montgomery_pow_kary.  It performs ARRAY_SIZE
// modular exponentiations when called.  It exists as an optimization that
// can have significantly higher throughput than the plain montgomery_pow_kary.
// It seems to be common to achieve the highest throughput with an ARRAY_SIZE of
// around 4, but of course you will need to benchmark to know what will be
// fastest for your particular machine and code.
//
template <class MF, typename U, size_t ARRAY_SIZE>
std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
array_montgomery_pow_kary(const std::array<MF, ARRAY_SIZE>& mf,
                          const std::array<typename MF::MontgomeryValue, ARRAY_SIZE>& x,
                          const std::array<U, ARRAY_SIZE>& n)
{
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    using V = typename MF::MontgomeryValue;
    using std::size_t;

    constexpr size_t P = 4;   // (1 << P) == the k in k-ary exponentiation

    // initialize the precalculation table for k-ary pow algorithm
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1 << P;   
    static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
    V table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        table[0][j] = mf[j].getUnityValue();   // montgomery one
        table[1][j] = x[j];
    }
    if (TABLESIZE >= 4) {
        // we should check for a ((power of 2) >= 4), but this is
        // probably adquate or our needs
        static_assert(TABLESIZE % 4 == 0, "");
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[2][j] = mf[j].square(x[j]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[3][j] = mf[j].multiply(x[j], table[2][j]);
        //HURCHALLA_REQUEST_UNROLL_LOOP
        for (std::size_t i=4; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i][j] = mf[j].square(table[halfi][j]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i+1][j] = mf[j].multiply(table[halfi][j], table[halfi+1][j]);
        }
    }

    U n_max = n[0];
    HURCHALLA_REQUEST_UNROLL_LOOP for (int j=1; j<ARRAY_SIZE; ++j) {
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
                result[j] = mf[j].square(result[j]);
        }
        shift -= static_cast<int>(P);
        std::array<size_t, ARRAY_SIZE> index;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
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
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mf[j].square(result[j]);
    }
    size_t tmpmask = (1 << shift) - 1;
    std::array<size_t, ARRAY_SIZE> index;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        index[j] = static_cast<size_t>(n[j]) & tmpmask;
        result[j] = mf[j].multiply(result[j], table[index[j]][j]);
    }
    return result;
}




// This is an alternative to calling MontgomeryForm's pow() member function.
// This function has a MontgomeryForm parameter (mf) since it is not a member
// function, but otherwise this function should be a drop-in replacement for
// the MontgomeryForm pow() member function.  For more information, see
// MontgomeryForm.h and the array version of the pow() function contained within
// it.
//
// Although I have not benchmarked this function yet, there is a very good
// chance that it will time faster than the pow() member function, for 128 bit
// types.
//
// Prospectively this function's capabilities will be entirely integrated into
// MontgomeryForm's pow() in the future, once this function is no longer
// experimental.
//
template <class MF, typename T, size_t ARRAY_SIZE>
std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
array_montgomery_pow_kary(const MF& mf,
                          const std::array<typename MF::MontgomeryValue, ARRAY_SIZE>& x,
                          T nexp)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    using U = typename extensible_make_unsigned<T>::type;

    HPBC_PRECONDITION(nexp >= 0);
    U n = static_cast<U>(nexp);

    using V = typename MF::MontgomeryValue;
    using std::size_t;

    constexpr size_t P = 4;   // (1 << P) == the k in k-ary exponentiation

    // initialize the precalculation table for k-ary pow algorithm
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1 << P;   
    static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
    V table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        table[0][j] = mf.getUnityValue();   // montgomery one
        table[1][j] = x[j];
    }
    if (TABLESIZE >= 4) {
        // we should check for a ((power of 2) >= 4), but this is
        // probably adquate or our needs
        static_assert(TABLESIZE % 4 == 0, "");
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[2][j] = mf.square(x[j]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[3][j] = mf.multiply(x[j], table[2][j]);
        //HURCHALLA_REQUEST_UNROLL_LOOP
        for (std::size_t i=4; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i][j] = mf.square(table[halfi][j]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i+1][j] = mf.multiply(table[halfi][j], table[halfi+1][j]);
        }
    }

    U n_max = n;

    constexpr size_t MASK = TABLESIZE - 1;
    if (n_max <= MASK) {
        std::array<V, ARRAY_SIZE> result;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = table[static_cast<size_t>(n)][j];
        return result;
    }


    // count_leading_zeros returns the number of leading 0-bits in n_max, starting
    // at the most significant bit position. If n_max is 0, the result is undefined
    HPBC_ASSERT(n_max > 0);
    int leading_zeros = count_leading_zeros(n_max);
    int numbits = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
    // because we returned above if (n_max <= MASK), we can assert the following:
    HPBC_ASSERT(numbits > static_cast<int>(P));

    int shift = numbits - static_cast<int>(P);
    std::array<V, ARRAY_SIZE> result;
    size_t tmp = static_cast<size_t>(n >> shift);
    HPBC_ASSERT(tmp <= MASK);
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
        result[j] = table[tmp][j];
    }


    while (shift >= static_cast<int>(P)) {
        for (size_t i=0; i<P; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf.square(result[j]);
        }
#if 1
// sliding window optimization
// TODO: maybe optimize the shift in the next line, since n may be > 64bit
        while (shift > static_cast<int>(P) && (static_cast<size_t>(n >> (shift-1)) & 1) == 0) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf.square(result[j]);
            --shift;
        }
#endif
        shift -= static_cast<int>(P);
// TODO: maybe optimize next line, since n may be > 64bit
        size_t index = static_cast<size_t>(n >> shift) & MASK;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf.multiply(result[j], table[index][j]);
        }
    }


    if (shift == 0)
        return result;
    HPBC_ASSERT(0 < shift && shift < static_cast<int>(P));

    for (int i=0; i<shift; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mf.square(result[j]);
    }
    size_t tmpmask = (1 << shift) - 1;
    size_t index = static_cast<size_t>(n) & tmpmask;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        result[j] = mf.multiply(result[j], table[index][j]);
    }
    return result;
}




} // end namespace

#endif
