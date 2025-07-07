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
// Implementation note: this is a modified version of the 2^k-ary exponentiation
// algorithm  ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring ),
// which precalculates the even exponents as well as the normal odd exponents,
// in order to avoid two conditional branches that would exist in the main loop
// of the normal 2^k-ary algorithm.  This is particularly helpful for the array
// version of this function further below.
//
template <class MF, typename T,
          bool USE_SLIDING_WINDOW_OPTIMIZATION = true,
          size_t TABLE_BITS = 4>
typename MF::MontgomeryValue montgomery_pow_kary(const MF& mf, typename MF::MontgomeryValue x, T nexp)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    // FYI you almost certainly want 2 <= TABLE_BITS <= 5.
    // TABLE_BITS > 0 is required; anything above 9 is probably a very bad idea
    // even if it works (it would mean calculating 1024+ table entries!)
    static_assert(0 < TABLE_BITS && TABLE_BITS < 10, "");

    HPBC_PRECONDITION(nexp >= 0);

    using U = typename extensible_make_unsigned<T>::type;
    using V = typename MF::MontgomeryValue;
    using std::size_t;
    U n = static_cast<U>(nexp);

    constexpr size_t P = TABLE_BITS;

    // initialize the precalculation table for 2^k-ary pow algorithm
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
        table[4] = mf.template square<LowuopsTag>(table[2]);
        table[5] = mf.template multiply<LowuopsTag>(table[2], table[3]);
        table[6] = mf.template square<LowuopsTag>(table[3]);
        table[7] = mf.template multiply<LowuopsTag>(table[3], table[4]);
    } else if (TABLESIZE == 16) {
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
        table[4] = mf.template square<LowuopsTag>(table[2]);
        table[5] = mf.template multiply<LowuopsTag>(table[2], table[3]);
        table[6] = mf.template square<LowuopsTag>(table[3]);
        table[7] = mf.template multiply<LowuopsTag>(table[3], table[4]);
        table[8] = mf.template square<LowuopsTag>(table[4]);
        table[9] = mf.template multiply<LowuopsTag>(table[4], table[5]);
        table[10] = mf.template square<LowuopsTag>(table[5]);
        table[11] = mf.template multiply<LowuopsTag>(table[5], table[6]);
        table[12] = mf.template square<LowuopsTag>(table[6]);
        table[13] = mf.template multiply<LowuopsTag>(table[6], table[7]);
        table[14] = mf.template square<LowuopsTag>(table[7]);
        table[15] = mf.template multiply<LowuopsTag>(table[7], table[8]);
    } else if (TABLESIZE == 32) {
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
        table[4] = mf.template square<LowuopsTag>(table[2]);
        table[5] = mf.template multiply<LowuopsTag>(table[2], table[3]);
        table[6] = mf.template square<LowuopsTag>(table[3]);
        table[7] = mf.template multiply<LowuopsTag>(table[3], table[4]);
        table[8] = mf.template square<LowuopsTag>(table[4]);
        table[9] = mf.template multiply<LowuopsTag>(table[4], table[5]);
        table[10] = mf.template square<LowuopsTag>(table[5]);
        table[11] = mf.template multiply<LowuopsTag>(table[5], table[6]);
        table[12] = mf.template square<LowuopsTag>(table[6]);
        table[13] = mf.template multiply<LowuopsTag>(table[6], table[7]);
        table[14] = mf.template square<LowuopsTag>(table[7]);
        table[15] = mf.template multiply<LowuopsTag>(table[7], table[8]);

        table[16] = mf.template square<LowuopsTag>(table[8]);
        table[17] = mf.template multiply<LowuopsTag>(table[8], table[9]);
        table[18] = mf.template square<LowuopsTag>(table[9]);
        table[19] = mf.template multiply<LowuopsTag>(table[9], table[10]);
        table[20] = mf.template square<LowuopsTag>(table[10]);
        table[21] = mf.template multiply<LowuopsTag>(table[10], table[11]);
        table[22] = mf.template square<LowuopsTag>(table[11]);
        table[23] = mf.template multiply<LowuopsTag>(table[11], table[12]);
        table[24] = mf.template square<LowuopsTag>(table[12]);
        table[25] = mf.template multiply<LowuopsTag>(table[12], table[13]);
        table[26] = mf.template square<LowuopsTag>(table[13]);
        table[27] = mf.template multiply<LowuopsTag>(table[13], table[14]);
        table[28] = mf.template square<LowuopsTag>(table[14]);
        table[29] = mf.template multiply<LowuopsTag>(table[14], table[15]);
        table[30] = mf.template square<LowuopsTag>(table[15]);
        table[31] = mf.template multiply<LowuopsTag>(table[15], table[16]);
    } else {
        // we should check for a ((power of 2) >= 64), but this is
        // probably adequate for our needs
        HPBC_ASSERT(TABLESIZE % 64 == 0);
        table[2] = mf.square(x);
        table[3] = mf.multiply(x, table[2]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=4; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            table[i] = mf.template square<LowuopsTag>(table[halfi]);
            table[i+1] = mf.template multiply<LowuopsTag>(table[halfi], table[halfi + 1]);
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






// This is the array version of montgomery_pow_kary.  It performs ARRAY_SIZE
// modular exponentiations when called.  It exists as an optimization that
// can have significantly higher throughput than the plain montgomery_pow_kary.
// It seems to be common to achieve the highest throughput with an ARRAY_SIZE of
// around 4, but of course you will need to benchmark to know what will be
// fastest for your particular machine and code.
//
template <class MF, typename U,
          size_t ARRAY_SIZE,
          size_t TABLE_BITS = 4>
std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
array_montgomery_pow_kary(const std::array<MF, ARRAY_SIZE>& mf,
                          const std::array<typename MF::MontgomeryValue, ARRAY_SIZE>& x,
                          const std::array<U, ARRAY_SIZE>& n)
{
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");
    // FYI you almost certainly want 2 <= TABLE_BITS <= 5.
    // TABLE_BITS > 0 is required; anything above 9 is probably a very bad idea
    // even if it works (it would mean calculating 1024+ table entries!)
    static_assert(0 < TABLE_BITS && TABLE_BITS < 10, "");

    using V = typename MF::MontgomeryValue;
    using std::size_t;

    constexpr size_t P = TABLE_BITS;

    // initialize the precalculation table for 2^k-ary pow algorithm
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1 << P;   
    static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
    V table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        table[0][j] = mf[j].getUnityValue();   // montgomery one
        table[1][j] = x[j];
    }
    if (TABLESIZE >= 4) {
        // we should ideally check for a ((power of 2) >= 4), but this assert
        // is probably good enough for our needs:
        HPBC_ASSERT(TABLESIZE % 4 == 0);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[2][j] = mf[j].template square<LowuopsTag>(x[j]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[3][j] = mf[j].template multiply<LowuopsTag>(x[j], table[2][j]);
        //HURCHALLA_REQUEST_UNROLL_LOOP
        for (std::size_t i=4; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i][j] = mf[j].template square<LowuopsTag>(table[halfi][j]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i+1][j] = mf[j].template multiply<LowuopsTag>(table[halfi][j], table[halfi+1][j]);
        }
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
template <class MF, typename T,
          size_t ARRAY_SIZE,
          bool USE_SLIDING_WINDOW_OPTIMIZATION = true,
          size_t TABLE_BITS = 4>
std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
array_montgomery_pow_kary(const MF& mf,
                          const std::array<typename MF::MontgomeryValue, ARRAY_SIZE>& x,
                          T nexp)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    // FYI you almost certainly want 2 <= TABLE_BITS <= 5.
    // TABLE_BITS > 0 is required; anything above 9 is probably a very bad idea
    // even if it works (it would mean calculating 1024+ table entries!)
    static_assert(0 < TABLE_BITS && TABLE_BITS < 10, "");

    HPBC_PRECONDITION(nexp >= 0);

    using U = typename extensible_make_unsigned<T>::type;
    using V = typename MF::MontgomeryValue;
    using std::size_t;
    U n = static_cast<U>(nexp);

    constexpr size_t P = TABLE_BITS;

    // initialize the precalculation table for 2^k-ary pow algorithm
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
            table[2][j] = mf.template square<LowuopsTag>(x[j]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[3][j] = mf.template multiply<LowuopsTag>(x[j], table[2][j]);
        //HURCHALLA_REQUEST_UNROLL_LOOP
        for (std::size_t i=4; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i][j] = mf.template square<LowuopsTag>(table[halfi][j]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i+1][j] = mf.template multiply<LowuopsTag>(table[halfi][j], table[halfi+1][j]);
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
        if (USE_SLIDING_WINDOW_OPTIMIZATION) {
            while (shift > static_cast<int>(P) && (static_cast<size_t>(n >> (shift-1)) & 1) == 0) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = mf.template square<LowuopsTag>(result[j]);
                --shift;
            }
        }
        for (size_t i=0; i<P; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf.template square<LowuopsTag>(result[j]);
        }
        shift -= static_cast<int>(P);
        size_t index = static_cast<size_t>(n >> shift) & MASK;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf.template multiply<LowuopsTag>(result[j], table[index][j]);
        }
    }


    if (shift == 0)
        return result;
    HPBC_ASSERT(0 < shift && shift < static_cast<int>(P));

    for (int i=0; i<shift; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mf.template square<LowuopsTag>(result[j]);
    }
    size_t tmpmask = (1 << shift) - 1;
    size_t index = static_cast<size_t>(n) & tmpmask;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        result[j] = mf.template multiply<LowuopsTag>(result[j], table[index][j]);
    }
    return result;
}




} // end namespace

#endif
