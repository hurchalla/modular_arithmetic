// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_MONTGOMERY_TWO_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_MONTGOMERY_TWO_POW_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontgomeryFormExtensions.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/count_leading_zeros.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/branchless_shift_right.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <type_traits>
#include <cstddef>
#include <array>
#include <utility>

namespace hurchalla { namespace detail {


// Implementation notes: this is a highly modified version of the 2^k-ary
// algorithm  ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring ),
// taking advantage of the fact that the base is always 2.
// It makes calls to twoPowLimited_times_x() to completely replace the
// ordinary table that would be used for 2^kary exponentiation, and it uses
// additional small (real) tables to further improve performance.
//
// All code here is copied from experimental/experimental_montgomery_two_pow.h.
// The CODE_SECTIONs used here are the exact same numbered CODE_SECTIONs in
// that file.  By reading initial comments in the corresponding CODE_SECTIONs
// inside that file, you can trace the origins of most CODE_SECTIONs back to
// referenced earlier CODE_SECTIONs in that file, which makes it possible to
// see changes and optimizations.
//
// Minor note: we use a struct with static member functions to disallow ADL.
struct impl_montgomery_two_pow {

  template <typename U>
  static constexpr int floor_log2(U x)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");
    // x > 0 required, but C++11 constexpr function won't allow this check
    //HPBC_CLOCKWORK_PRECONDITION2(x > 0);
    return (x <= 1) ? 0 : 1 + impl_montgomery_two_pow::floor_log2(x >> 1);
  }


  // Calculate pow(2, n), modulo the modulus of mf, and return the result in
  // montgomeryform representation.
  template <class MF, typename U,
            bool USE_SLIDING_WINDOW_OPTIMIZATION,
            size_t TABLE_BITS, size_t CODE_SECTION,
            bool USE_SQUARING_VALUE_OPTIMIZATION>
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // note: The unused template parameter TABLE_BITS exists simply to maintain
    // correspondence with experimental_montgomery_two_pow.h.
    //
    // All of this code is a copy of experimental_montgomery_two_pow.h's
    // non-array call(), for TABLE_BITS == 0 with the following code sections.
    static_assert(TABLE_BITS == 0 &&
           ((CODE_SECTION >= 22 && CODE_SECTION <= 26) ||
            (CODE_SECTION >= 31 && CODE_SECTION <= 33) ||
            (CODE_SECTION >= 34 && CODE_SECTION <= 38)), "");

    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    namespace hc = hurchalla;
    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using MFE = hc::detail::MontgomeryFormExtensions<MF, hc::LowlatencyTag>;
    using SV = typename MFE::SquaringValue;
    using std::size_t;

    using RU = typename MFE::RU;
    constexpr int digitsRU = hc::ut_numeric_limits<RU>::digits;
    constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));
    constexpr size_t MASK = (1u << P2) - 1u;

    using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;


    if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 22 && CODE_SECTION <= 26) {
        // this section corresponds exactly to CODE_SECTIONs 22 - 26 in
        // experimental_montgomery_two_pow.h.  It should be essentially a
        // copy/paste.

        constexpr int NUMBITS_TABLE_HIGH_SIZE = 2;
        constexpr size_t NUM_EXTRA_TABLES = 2 * (CODE_SECTION - 21);

        constexpr size_t TABLE_HIGH_SIZE = 1u << NUMBITS_TABLE_HIGH_SIZE;
        constexpr int P3 = P2 + NUMBITS_TABLE_HIGH_SIZE;
        constexpr int NUMBITS_MASKBIG = P3 + NUM_EXTRA_TABLES * NUMBITS_TABLE_HIGH_SIZE;
        static_assert(std::numeric_limits<size_t>::digits > NUMBITS_MASKBIG, "");
        constexpr size_t MASKBIG = (static_cast<size_t>(1) << NUMBITS_MASKBIG) - 1u;
        std::array<C, TABLE_HIGH_SIZE> table_mid;
        std::array<std::array<V, TABLE_HIGH_SIZE>, NUM_EXTRA_TABLES> tables_extra;

        int shift = 0;
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
        }

        V result;
        static_assert(TABLE_HIGH_SIZE == 4, "");

        {
            C r1 = MFE::getMontvalueR(mf);
            V r2 = mf.square(r1);
            V r3 = mf.template multiply<LowuopsTag>(r2, r1);
            V r4 = mf.square(r2);
            table_mid[0] = r1;                         // R^1
            table_mid[1] = mf.getCanonicalValue(r2);   // R^2
            table_mid[2] = mf.getCanonicalValue(r3);   // R^3
            table_mid[3] = mf.getCanonicalValue(r4);   // R^4

            HPBC_CLOCKWORK_ASSERT2(shift >= 0);
            size_t tmp = static_cast<size_t>(branchless_shift_right(n, shift));
            HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
            size_t loindex = tmp & MASK;
            size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            result = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            V next = r4;              // R^4
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                tables_extra[i][0] = mf.getUnityValue();   // R^0
                tables_extra[i][1] = next;
                V nextSq = mf.square(next);
                V nexttmp = mf.square(nextSq);
                tables_extra[i][2] = nextSq;
                tables_extra[i][3] = mf.template multiply<LowuopsTag>(nextSq, next);
                next = nexttmp;

                int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                result = mf.template multiply<LowuopsTag>(tables_extra[i][index_extra], result);
            }
        }


        while (shift >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (shift > NUMBITS_MASKBIG && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                        sv = MFE::squareSV(mf, sv);
                        --shift;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

                shift -= NUMBITS_MASKBIG;
                size_t tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                size_t loindex = tmp & MASK;
                size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
                V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P3-1; ++i)
                    sv = MFE::squareSV(mf, sv);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                    int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                    size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                    val1 = mf.template multiply<LowuopsTag>(val1, tables_extra[i][index_extra]);

                    static_assert(NUMBITS_TABLE_HIGH_SIZE == 2, "");
                    sv = MFE::squareSV(mf, sv);
                    sv = MFE::squareSV(mf, sv);
                }
                result = MFE::squareToMontgomeryValue(mf, sv);

                result = mf.multiply(result, val1);
            }
            else {
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (shift > NUMBITS_MASKBIG && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                        result = mf.square(result);
                        --shift;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

                shift -= NUMBITS_MASKBIG;
                size_t tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                size_t loindex = tmp & MASK;
                size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
                V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P3; ++i)
                    result = mf.square(result);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                    int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                    size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                    val1 = mf.template multiply<LowuopsTag>(val1, tables_extra[i][index_extra]);

                    static_assert(NUMBITS_TABLE_HIGH_SIZE == 2, "");
                    result = mf.square(result);
                    result = mf.square(result);
                }

                result = mf.multiply(result, val1);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        size_t tmp = static_cast<size_t>(n) & tmpmask;
        size_t loindex = tmp & MASK;
        size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
        V val1 = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            // could use:
            //for (int i=0; i < NUM_EXTRA_TABLES && (2*i + P3 < shift); ++i)

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                val1 = mf.multiply(val1, tables_extra[i][index_extra]);
            }

            SV sv = MFE::getSquaringValue(mf, result);
            HPBC_CLOCKWORK_ASSERT2(shift >= 1);
            for (int i=0; i<shift-1; ++i)
                sv = MFE::squareSV(mf, sv);
            result = MFE::squareToMontgomeryValue(mf, sv);
        }
        else {
            result = mf.square(result);

            // could use:
            //for (int i=0; i < NUM_EXTRA_TABLES && (2*i + P3 < shift); ++i)

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                val1 = mf.multiply(val1, tables_extra[i][index_extra]);
            }

            for (int i=1; i<shift; ++i)
                result = mf.square(result);
        }
        result = mf.multiply(result, val1);
        return result;

    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 31 && CODE_SECTION <= 33) {
        // this section corresponds exactly to CODE_SECTIONs 31 - 33 in
        // experimental_montgomery_two_pow.h.  It should be essentially a
        // copy/paste.

        // make the high tables' size either 2, 4, 8, depending on the CODE_SECTION
        constexpr int NUMBITS_TABLE_HIGH_SIZE = static_cast<int>(CODE_SECTION) - 30;

        static_assert(1 <= NUMBITS_TABLE_HIGH_SIZE && NUMBITS_TABLE_HIGH_SIZE <= 3, "");

        constexpr size_t TABLE_HIGH_SIZE = 1u << NUMBITS_TABLE_HIGH_SIZE;
        constexpr int NUMBITS_MASKBIG = P2 + 3 * NUMBITS_TABLE_HIGH_SIZE;
        constexpr int P3 = P2 + NUMBITS_TABLE_HIGH_SIZE;
        constexpr int P4 = P3 + NUMBITS_TABLE_HIGH_SIZE;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;
        std::array<C, TABLE_HIGH_SIZE> table_mid;
        std::array<V, TABLE_HIGH_SIZE> table_high;
        std::array<V, TABLE_HIGH_SIZE> table_ultra;

        int shift = 0;
        size_t tmp;
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
            HPBC_CLOCKWORK_ASSERT2(shift >= 0);
            tmp = static_cast<size_t>(branchless_shift_right(n, shift));
            n = branchless_shift_left(n, leading_zeros + NUMBITS_MASKBIG);
        }
        else {
            HPBC_CLOCKWORK_ASSERT2(n <= MASKBIG);
            tmp = static_cast<size_t>(n);
        }

        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        size_t loindex = tmp & MASK;
        size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
        size_t hiindex = (tmp >> P3) & (TABLE_HIGH_SIZE - 1);
        size_t ultindex = (tmp >> P4);
        HPBC_CLOCKWORK_ASSERT2(ultindex < TABLE_HIGH_SIZE);

        V result;

        static_assert(2 <= TABLE_HIGH_SIZE && TABLE_HIGH_SIZE <= 8, "");
        if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 2) {
            C r1 = MFE::getMontvalueR(mf);
            V r2 = mf.square(r1);
            table_mid[0] = r1;                             // R^1
            table_mid[1] = mf.getCanonicalValue(r2);       // R^2

            table_high[0] = mf.getUnityValue();            // R^0
            table_high[1] = r2;                            // R^2

            V ttx = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            table_ultra[0] = mf.getUnityValue();           // R^0
            table_ultra[1] = mf.square(r2);                // R^4

            V val1 = mf.multiply(ttx, table_high[hiindex]);
            result = mf.multiply(val1, table_ultra[ultindex]);
        }
        else if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 4) {
            C r1 = MFE::getMontvalueR(mf);
            V r2 = mf.square(r1);
            V r3 = mf.template multiply<LowuopsTag>(r2, r1);
            V r4 = mf.square(r2);
            table_mid[0] = r1;                         // R^1
            table_mid[1] = mf.getCanonicalValue(r2);   // R^2
            table_mid[2] = mf.getCanonicalValue(r3);   // R^3
            table_mid[3] = mf.getCanonicalValue(r4);   // R^4

            V ttx = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            V r8 = mf.square(r4);
            V r16 = mf.square(r8);

            table_high[0] = mf.getUnityValue();      // R^0
            table_high[1] = r4;                      // R^4
            table_high[2] = r8;                      // R^8
            table_high[3] = mf.multiply(r8, r4);     // R^12

            V r32 = mf.square(r16);

            V val1 = mf.multiply(ttx, table_high[hiindex]);

            table_ultra[0] = mf.getUnityValue();     // R^0
            table_ultra[1] = r16;                    // R^16
            table_ultra[2] = r32;                    // R^32
            table_ultra[3] = mf.multiply(r32, r16);  // R^48

            result = mf.multiply(val1, table_ultra[ultindex]);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 8) {
            C r1 = MFE::getMontvalueR(mf);
            V r2 = mf.square(r1);
            V r4 = mf.square(r2);
            V r3 = mf.multiply(r2, r1);

            table_mid[0] = r1;                                                             // R^1
            table_mid[1] = mf.getCanonicalValue(r2);                                       // R^2

            V r8 = mf.square(r4);

            table_mid[2] = mf.getCanonicalValue(r3);                                       // R^3
            table_mid[3] = mf.getCanonicalValue(r4);                                       // R^4
            table_mid[4] = mf.getCanonicalValue(mf.template multiply<LowuopsTag>(r3, r2)); // R^5

            V r16 = mf.square(r8);

            table_mid[5] = mf.getCanonicalValue(mf.template square<LowuopsTag>(r3));       // R^6
            table_mid[6] = mf.getCanonicalValue(mf.template multiply<LowuopsTag>(r3, r4)); // R^7
            table_mid[7] = mf.getCanonicalValue(r8);                                       // R^8

            V r32 = mf.square(r16);
            V r24 = mf.multiply(r16, r8);

            V ttx = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            table_high[0] = mf.getUnityValue();                          // R^0
            table_high[1] = r8;                                          // R^8
            table_high[2] = r16;                                         // R^16
            table_high[3] = r24;                                         // R^24

            V r64 = mf.square(r32);

            table_high[4] = r32;                                         // R^32
            table_high[5] = mf.template multiply<LowuopsTag>(r24, r16);  // R^40

            V r128 = mf.square(r64);
            table_high[6] = mf.template square<LowuopsTag>(r24);         // R^48
            table_high[7] = mf.template multiply<LowuopsTag>(r32, r24);  // R^56

            V r192 = mf.multiply(r128, r64);
            V r256 = mf.square(r128);

            V val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

            table_ultra[0] = mf.getUnityValue();                            // R^0
            table_ultra[1] = r64;                                           // R^64
            table_ultra[2] = r128;                                          // R^128
            table_ultra[3] = r192;                                          // R^192
            table_ultra[4] = r256;                                          // R^256
            table_ultra[5] = mf.template multiply<LowuopsTag>(r192, r128);  // R^320
            table_ultra[6] = mf.square(r192);                               // R^384
            table_ultra[7] = mf.multiply(r256, r192);                       // R^448

            result = mf.multiply(val1, table_ultra[ultindex]);
        }

        // note: 'n' is type U.
        // Calculate the constexpr var 'high_word_shift' - when we right shift a
        // type U variable by this amount, we'll get the size_t furthest most
        // left bits of the type U variable.  Note that we assume that a right
        // shift by high_word_shift will be zero cost, since the shift is just a
        // way to access the CPU register that has the most significant bits -
        // unless the compiler is really dumb and misses this optimization,
        // which I haven't seen happen and which would surprise me.
        constexpr int size_t_digits = ut_numeric_limits<size_t>::digits;
        constexpr int digits_U = ut_numeric_limits<U>::digits;
        constexpr int digits_bigger = (digits_U > size_t_digits) ? digits_U : size_t_digits;
        constexpr int digits_smaller = (digits_U < size_t_digits) ? digits_U : size_t_digits;
        constexpr int high_word_shift = digits_bigger - size_t_digits;
        // the conditional below is just to avoid a compiler warning about a
        // negative shift in the loop, even though it would never happen
        constexpr int small_shift = (digits_smaller < NUMBITS_MASKBIG)
                                     ? 0 : (digits_smaller - NUMBITS_MASKBIG);

        int bits_remaining = shift;

        while (bits_remaining >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (bits_remaining > NUMBITS_MASKBIG &&
                                (static_cast<size_t>(n >> high_word_shift) &
                                     (static_cast<size_t>(1) << (digits_smaller - 1))) == 0) {
                        sv = MFE::squareSV(mf, sv);
                        n = static_cast<U>(n << 1);
                        --bits_remaining;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(bits_remaining >= NUMBITS_MASKBIG);

//                tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
                n = static_cast<U>(n << NUMBITS_MASKBIG);
                bits_remaining -= NUMBITS_MASKBIG;

                loindex = tmp & MASK;
                midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
                hiindex = (tmp >> P3) & (TABLE_HIGH_SIZE - 1);
                ultindex = (tmp >> P4) & (TABLE_HIGH_SIZE - 1);

                V ttx = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

                static_assert(NUMBITS_MASKBIG >= 5, "");
                sv = MFE::squareSV(mf, sv);
                sv = MFE::squareSV(mf, sv);

                V val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

                sv = MFE::squareSV(mf, sv);
                sv = MFE::squareSV(mf, sv);

                val1 = mf.template multiply<LowuopsTag>(val1, table_ultra[ultindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=4; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);

                result = mf.multiply(result, val1);
            }
            else {
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (bits_remaining > NUMBITS_MASKBIG &&
                                (static_cast<size_t>(n >> high_word_shift) &
                                     (static_cast<size_t>(1) << (digits_smaller - 1))) == 0) {
                        result = mf.square(result);
                        n = static_cast<U>(n << 1);
                        --bits_remaining;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(bits_remaining >= NUMBITS_MASKBIG);

                //tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
                n = static_cast<U>(n << NUMBITS_MASKBIG);
                bits_remaining -= NUMBITS_MASKBIG;

                loindex = tmp & MASK;
                midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
                hiindex = (tmp >> P3) & (TABLE_HIGH_SIZE - 1);
                ultindex = (tmp >> P4) & (TABLE_HIGH_SIZE - 1);

                V ttx = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

                static_assert(NUMBITS_MASKBIG >= 4, "");
                result = mf.square(result);
                result = mf.square(result);

                V val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

                result = mf.square(result);
                result = mf.square(result);

                val1 = mf.template multiply<LowuopsTag>(val1, table_ultra[ultindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=4; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);

                result = mf.multiply(result, val1);
            }
        }
        if (bits_remaining == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < bits_remaining && bits_remaining < NUMBITS_MASKBIG);

        //size_t tmpmask = (1u << shift) - 1u;
        //size_t index = static_cast<size_t>(n) & tmpmask;
        size_t index = static_cast<size_t>(n >> high_word_shift) >> (digits_smaller - bits_remaining);
        loindex = index & MASK;
        midindex = (index >> P2) & (TABLE_HIGH_SIZE - 1);
        hiindex = (index >> P3) & (TABLE_HIGH_SIZE - 1);
        ultindex = (index >> P4) & (TABLE_HIGH_SIZE - 1);
        V ttx = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

        V val1;
        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);
            val1 = mf.template multiply<LowuopsTag>(val1, table_ultra[ultindex]);

            SV sv = MFE::getSquaringValue(mf, result);
            HPBC_CLOCKWORK_ASSERT2(bits_remaining >= 1);
            for (int i=0; i<bits_remaining-1; ++i)
                sv = MFE::squareSV(mf, sv);
            result = MFE::squareToMontgomeryValue(mf, sv);
        }
        else {
            result = mf.square(result);
            val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);
            val1 = mf.template multiply<LowuopsTag>(val1, table_ultra[ultindex]);

            for (int i=1; i<bits_remaining; ++i)
                result = mf.square(result);
        }
        result = mf.multiply(result, val1);
        return result;

    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 34 && CODE_SECTION <= 38) {
        // this section corresponds exactly to CODE_SECTIONs 34 - 38 in
        // experimental_montgomery_two_pow.h.  It should be essentially a
        // copy/paste.

        constexpr int NUMBITS_TABLE_HIGH_SIZE = 2;
        constexpr size_t NUM_EXTRA_TABLES = 2 * (CODE_SECTION - 33);

        constexpr size_t TABLE_HIGH_SIZE = 1u << NUMBITS_TABLE_HIGH_SIZE;
        constexpr int P3 = P2 + NUMBITS_TABLE_HIGH_SIZE;
        constexpr int NUMBITS_MASKBIG = P3 + NUM_EXTRA_TABLES * NUMBITS_TABLE_HIGH_SIZE;
        static_assert(std::numeric_limits<size_t>::digits > NUMBITS_MASKBIG, "");
        constexpr size_t MASKBIG = (static_cast<size_t>(1) << NUMBITS_MASKBIG) - 1u;
        std::array<C, TABLE_HIGH_SIZE> table_mid;
        std::array<std::array<V, TABLE_HIGH_SIZE>, NUM_EXTRA_TABLES> tables_extra;

        int shift = 0;
        size_t tmp;
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
            HPBC_CLOCKWORK_ASSERT2(shift >= 0);
            tmp = static_cast<size_t>(branchless_shift_right(n, shift));
            n = branchless_shift_left(n, leading_zeros + NUMBITS_MASKBIG);
        }
        else {
            HPBC_CLOCKWORK_ASSERT2(n <= MASKBIG);
            tmp = static_cast<size_t>(n);
        }

        V result;
        static_assert(TABLE_HIGH_SIZE == 4, "");

        {
            C r1 = MFE::getMontvalueR(mf);
            V r2 = mf.square(r1);
            V r3 = mf.template multiply<LowuopsTag>(r2, r1);
            V r4 = mf.square(r2);
            table_mid[0] = r1;                         // R^1
            table_mid[1] = mf.getCanonicalValue(r2);   // R^2
            table_mid[2] = mf.getCanonicalValue(r3);   // R^3
            table_mid[3] = mf.getCanonicalValue(r4);   // R^4

            HPBC_CLOCKWORK_ASSERT2(shift >= 0);
            //size_t tmp = static_cast<size_t>(branchless_shift_right(n, shift));
            HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
            size_t loindex = tmp & MASK;
            size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            result = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            V next = r4;              // R^4
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                tables_extra[i][0] = mf.getUnityValue();   // R^0
                tables_extra[i][1] = next;
                V nextSq = mf.square(next);
                V nexttmp = mf.square(nextSq);
                tables_extra[i][2] = nextSq;
                tables_extra[i][3] = mf.template multiply<LowuopsTag>(nextSq, next);
                next = nexttmp;

                int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                result = mf.template multiply<LowuopsTag>(tables_extra[i][index_extra], result);
            }
        }

        // note: 'n' is type U.
        // Calculate the constexpr var 'high_word_shift' - when we right shift a
        // type U variable by this amount, we'll get the size_t furthest most
        // left bits of the type U variable.  Note that we assume that a right
        // shift by high_word_shift will be zero cost, since the shift is just a
        // way to access the CPU register that has the most significant bits -
        // unless the compiler is really dumb and misses this optimization,
        // which I haven't seen happen and which would surprise me.
        constexpr int size_t_digits = ut_numeric_limits<size_t>::digits;
        constexpr int digits_U = ut_numeric_limits<U>::digits;
        constexpr int digits_bigger = (digits_U > size_t_digits) ? digits_U : size_t_digits;
        constexpr int digits_smaller = (digits_U < size_t_digits) ? digits_U : size_t_digits;
        constexpr int high_word_shift = digits_bigger - size_t_digits;
        // the conditional below is just to avoid a compiler warning about a
        // negative shift in the loop, even though it would never happen
        constexpr int small_shift = (digits_smaller < NUMBITS_MASKBIG)
                                     ? 0 : (digits_smaller - NUMBITS_MASKBIG);

        int bits_remaining = shift;


        while (bits_remaining >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (bits_remaining > NUMBITS_MASKBIG &&
                                (static_cast<size_t>(n >> high_word_shift) &
                                     (static_cast<size_t>(1) << (digits_smaller - 1))) == 0) {
                        sv = MFE::squareSV(mf, sv);
                        n = static_cast<U>(n << 1);
                        --bits_remaining;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(bits_remaining >= NUMBITS_MASKBIG);

                //tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
                n = static_cast<U>(n << NUMBITS_MASKBIG);
                bits_remaining -= NUMBITS_MASKBIG;

                size_t loindex = tmp & MASK;
                size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
                V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P3-1; ++i)
                    sv = MFE::squareSV(mf, sv);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                    int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                    size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                    val1 = mf.template multiply<LowuopsTag>(val1, tables_extra[i][index_extra]);

                    static_assert(NUMBITS_TABLE_HIGH_SIZE == 2, "");
                    sv = MFE::squareSV(mf, sv);
                    sv = MFE::squareSV(mf, sv);
                }
                result = MFE::squareToMontgomeryValue(mf, sv);

                result = mf.multiply(result, val1);
            }
            else {
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (bits_remaining > NUMBITS_MASKBIG &&
                                (static_cast<size_t>(n >> high_word_shift) &
                                     (static_cast<size_t>(1) << (digits_smaller - 1))) == 0) {
                        result = mf.square(result);
                        n = static_cast<U>(n << 1);
                        --bits_remaining;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(bits_remaining >= NUMBITS_MASKBIG);

                //tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
                n = static_cast<U>(n << NUMBITS_MASKBIG);
                bits_remaining -= NUMBITS_MASKBIG;

                size_t loindex = tmp & MASK;
                size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
                V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P3; ++i)
                    result = mf.square(result);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                    int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                    size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                    val1 = mf.template multiply<LowuopsTag>(val1, tables_extra[i][index_extra]);

                    static_assert(NUMBITS_TABLE_HIGH_SIZE == 2, "");
                    result = mf.square(result);
                    result = mf.square(result);
                }

                result = mf.multiply(result, val1);
            }
        }
        if (bits_remaining == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < bits_remaining && bits_remaining < NUMBITS_MASKBIG);

        //size_t tmpmask = (1u << shift) - 1u;
        //tmp = static_cast<size_t>(n) & tmpmask;
        tmp = static_cast<size_t>(n >> high_word_shift) >> (digits_smaller - bits_remaining);
        size_t loindex = tmp & MASK;
        size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
        V val1 = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            // could use:
            //for (int i=0; i < NUM_EXTRA_TABLES && (2*i + P3 < bits_remaining); ++i)

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                val1 = mf.multiply(val1, tables_extra[i][index_extra]);
            }

            SV sv = MFE::getSquaringValue(mf, result);
            HPBC_CLOCKWORK_ASSERT2(bits_remaining >= 1);
            for (int i=0; i<bits_remaining-1; ++i)
                sv = MFE::squareSV(mf, sv);
            result = MFE::squareToMontgomeryValue(mf, sv);
        }
        else {
            result = mf.square(result);

            // could use:
            //for (int i=0; i < NUM_EXTRA_TABLES && (2*i + P3 < bits_remaining); ++i)

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i < NUM_EXTRA_TABLES; ++i) {
                int P_extra = P3 + static_cast<int>(i * NUMBITS_TABLE_HIGH_SIZE);
                size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                val1 = mf.multiply(val1, tables_extra[i][index_extra]);
            }

            for (int i=1; i<bits_remaining; ++i)
                result = mf.square(result);
        }
        result = mf.multiply(result, val1);
        return result;
    }

  }  // end of function






#ifdef __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpass-failed"
#endif


  // Array version of montgomery two pow
  template <class MF, typename U,
            size_t ARRAY_SIZE, size_t TABLE_BITS, size_t CODE_SECTION,
            bool USE_SQUARING_VALUE_OPTIMIZATION>
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    namespace hc = hurchalla;
    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;
    using SV = typename MFE_LU::SquaringValue;
    using std::size_t;

    using RU = typename MFE_LU::RU;
    constexpr int digitsRU = hc::ut_numeric_limits<RU>::digits;
    constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));
    constexpr size_t MASK = (1u << P2) - 1u;

    U n_max = n[0];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=1; j<ARRAY_SIZE; ++j) {
        n_max = (n_max < n[j]) ? n[j] : n_max;
    }

    // note: The unused template parameter TABLE_BITS exists simply to maintain
    // correspondence with experimental_montgomery_two_pow.cpp.

    // note: This code is a copy of experimental_montgomery_two_pow.cpp's
    // array call(), for TABLE_BITS == 0, CODE_SECTIONs 28 through 31.

    HPBC_CLOCKWORK_ASSERT2(TABLE_BITS == 0 &&
                           (CODE_SECTION >= 28 && CODE_SECTION <= 31));

    if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 28) {
        std::array<V, ARRAY_SIZE> result;
        if (n_max <= MASK) {
            for (size_t j=0; j<ARRAY_SIZE; ++j) {
                C cR1 = MFE_LU::getMontvalueR(mf[j]);
                result[j]= MFE_LU::twoPowLimited_times_x(mf[j], static_cast<size_t>(n[j]), cR1);
            }
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);

        int shift = numbits - P2;
        HPBC_CLOCKWORK_ASSERT2(shift > 0);

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(branchless_shift_right(n[j], shift));
            HPBC_CLOCKWORK_ASSERT2(index <= MASK);
            // normally we use (index & MASK), but it's redundant with index <= MASK
            C cR1 = MFE_LU::getMontvalueR(mf[j]);
            result[j] = MFE_LU::twoPowLimited_times_x_v2(mf[j], index + 1, cR1);
        }

        while (shift >= P2) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                std::array<SV, ARRAY_SIZE> sv;
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::getSquaringValue(mf[j], result[j]);
                static_assert(P2 > 0, "");
                for (int i=0; i<P2 - 1; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf[j], sv[j]);
                }
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = MFE_LU::squareToMontgomeryValue(mf[j], sv[j]);
            } else {
                for (int i=0; i<P2; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                }
            }

            shift -= P2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t index = static_cast<size_t>(branchless_shift_right(n[j], shift)) & MASK;
                C tmp = mf[j].getCanonicalValue(result[j]);
                result[j] = MFE_LU::twoPowLimited_times_x_v2(mf[j], index + 1, tmp);
            }
        }

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].halve(result[j]);
        }

        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);
        size_t tmpmask = (1u << shift) - 1u;

        std::array<V, ARRAY_SIZE> tableVal;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n[j]) & tmpmask;
            C cR1 = MFE_LU::getMontvalueR(mf[j]);
            tableVal[j] = MFE_LU::twoPowLimited_times_x(mf[j], index, cR1);
        }
        for (int i=0; i<shift; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   tableVal[j]);
        }
        return result;
    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 29) {
        std::array<V, ARRAY_SIZE> result;
        if (n_max <= MASK) {
            for (size_t j=0; j<ARRAY_SIZE; ++j) {
                C cR1 = MFE_LU::getMontvalueR(mf[j]);
                result[j]= MFE_LU::twoPowLimited_times_x(mf[j], static_cast<size_t>(n[j]), cR1);
            }
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);

        int shift = numbits - P2;
        HPBC_CLOCKWORK_ASSERT2(shift > 0);

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = MFE_LU::getMontvalueR(mf[j]);
        }

        while (shift >= P2) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t index = static_cast<size_t>(branchless_shift_right(n[j], shift)) & MASK;
                C tmp = mf[j].getCanonicalValue(result[j]);
                result[j] = MFE_LU::twoPowLimited_times_x_v2(mf[j], index + 1, tmp);
            }

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                std::array<SV, ARRAY_SIZE> sv;
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::getSquaringValue(mf[j], result[j]);
                static_assert(P2 > 0, "");
                for (int i=0; i<P2 - 1; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf[j], sv[j]);
                }
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = MFE_LU::squareToMontgomeryValue(mf[j], sv[j]);
            } else {
                for (int i=0; i<P2; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                }
            }

            shift -= P2;
        }

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(branchless_shift_right(n[j], shift)) & MASK;
            C tmp = mf[j].getCanonicalValue(result[j]);
            result[j] = MFE_LU::twoPowLimited_times_x(mf[j], index, tmp);
        }

        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);
        size_t tmpmask = (1u << shift) - 1u;

        std::array<V, ARRAY_SIZE> tableVal;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n[j]) & tmpmask;
            C cR1 = MFE_LU::getMontvalueR(mf[j]);
            tableVal[j] = MFE_LU::twoPowLimited_times_x(mf[j], index, cR1);
        }
        for (int i=0; i<shift; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   tableVal[j]);
        }
        return result;
    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 30) {
        std::array<V, ARRAY_SIZE> result;
        if (n_max <= MASK) {
            for (size_t j=0; j<ARRAY_SIZE; ++j) {
                C cR1 = MFE_LU::getMontvalueR(mf[j]);
                result[j]= MFE_LU::twoPowLimited_times_x(mf[j], static_cast<size_t>(n[j]), cR1);
            }
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int bits_remaining = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(bits_remaining > P2);

        std::array<U, ARRAY_SIZE> n2;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
              // n2[j] = n[j] << leading_zeros;
            n2[j] = branchless_shift_left(n[j], leading_zeros);
        }

        // calculate the constexpr var 'high_word_shift' - when we right shift a
        // type U variable by this amount, we'll get the size_t furthest most
        // left bits of the type U variable.  Note that we assume that a right
        // shift by high_word_shift will be zero cost, since the shift is just a
        // way to access the CPU register that has the most significant bits -
        // unless the compiler is really dumb and misses this optimization,
        // which I haven't seen happen and which would surprise me.
        constexpr int size_t_digits = ut_numeric_limits<size_t>::digits;
        constexpr int digits_U = ut_numeric_limits<U>::digits;
        constexpr int digits_bigger = (digits_U > size_t_digits) ? digits_U : size_t_digits;
        constexpr int digits_smaller = (digits_U < size_t_digits) ? digits_U : size_t_digits;
        constexpr int high_word_shift = digits_bigger - size_t_digits;

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller - P2);
            n2[j] = static_cast<U>(n2[j] << P2);
            HPBC_CLOCKWORK_ASSERT2(index <= MASK);
            // normally we use (index & MASK), but it's redundant with index <= MASK
            C cR1 = MFE_LU::getMontvalueR(mf[j]);
            result[j] = MFE_LU::twoPowLimited_times_x_v2(mf[j], index + 1, cR1);
        }
        bits_remaining -= P2;

        while (bits_remaining >= P2) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                std::array<SV, ARRAY_SIZE> sv;
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::getSquaringValue(mf[j], result[j]);
                static_assert(P2 > 0, "");
                for (int i=0; i<P2 - 1; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf[j], sv[j]);
                }
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = MFE_LU::squareToMontgomeryValue(mf[j], sv[j]);
            } else {
                for (int i=0; i<P2; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                }
            }

            bits_remaining -= P2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller - P2);
                n2[j] = static_cast<U>(n2[j] << P2);
                C tmp = mf[j].getCanonicalValue(result[j]);
                result[j] = MFE_LU::twoPowLimited_times_x_v2(mf[j], index + 1, tmp);
            }
        }

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].halve(result[j]);
        }

        if (bits_remaining == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < bits_remaining && bits_remaining < P2);

        std::array<V, ARRAY_SIZE> tableVal;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller - bits_remaining);
            C cR1 = MFE_LU::getMontvalueR(mf[j]);
            tableVal[j] = MFE_LU::twoPowLimited_times_x(mf[j], index, cR1);
        }
        for (int i=0; i<bits_remaining; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   tableVal[j]);
        }
        return result;
    } else {      // CODE_SECTION 31
        static_assert(CODE_SECTION == 31, "");

        std::array<V, ARRAY_SIZE> result;
        if (n_max <= MASK) {
            for (size_t j=0; j<ARRAY_SIZE; ++j) {
                C cR1 = MFE_LU::getMontvalueR(mf[j]);
                result[j]= MFE_LU::twoPowLimited_times_x(mf[j], static_cast<size_t>(n[j]), cR1);
            }
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int bits_remaining = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(bits_remaining > P2);

        std::array<U, ARRAY_SIZE> n2;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
              // n2[j] = n[j] << leading_zeros;
            n2[j] = branchless_shift_left(n[j], leading_zeros);
        }

        // calculate the constexpr var 'high_word_shift' - when we right shift a
        // type U variable by this amount, we'll get the size_t furthest most
        // left bits of the type U variable.  Note that we assume that a right
        // shift by high_word_shift will be zero cost, since the shift is just a
        // way to access the CPU register that has the most significant bits -
        // unless the compiler is really dumb and misses this optimization,
        // which I haven't seen happen and which would surprise me.
        constexpr int size_t_digits = ut_numeric_limits<size_t>::digits;
        constexpr int digits_U = ut_numeric_limits<U>::digits;
        constexpr int digits_bigger = (digits_U > size_t_digits) ? digits_U : size_t_digits;
        constexpr int digits_smaller = (digits_U < size_t_digits) ? digits_U : size_t_digits;
        constexpr int high_word_shift = digits_bigger - size_t_digits;

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = MFE_LU::getMontvalueR(mf[j]);
        }

        HPBC_CLOCKWORK_ASSERT2(bits_remaining > P2);
        // we check against P2 + P2 because we always process P2 more bits after
        // the loop ends -- so we need to ensure we'll actually have
        // (bits_remaining >= P2) after the loop ends.
        while (bits_remaining >= P2 + P2) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller - P2);
                n2[j] = static_cast<U>(n2[j] << P2);
                C tmp = mf[j].getCanonicalValue(result[j]);
                result[j] = MFE_LU::twoPowLimited_times_x_v2(mf[j], index + 1, tmp);
            }

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                std::array<SV, ARRAY_SIZE> sv;
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::getSquaringValue(mf[j], result[j]);
                static_assert(P2 > 0, "");
                for (int i=0; i<P2 - 1; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf[j], sv[j]);
                }
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = MFE_LU::squareToMontgomeryValue(mf[j], sv[j]);
            } else {
                for (int i=0; i<P2; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                }
            }

            bits_remaining -= P2;
        }
        HPBC_CLOCKWORK_ASSERT2(P2 <= bits_remaining && bits_remaining < P2 + P2);

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller - P2);
            n2[j] = static_cast<U>(n2[j] << P2);

            C tmp = mf[j].getCanonicalValue(result[j]);
            result[j] = MFE_LU::twoPowLimited_times_x(mf[j], index, tmp);
        }
        bits_remaining -= P2;

        if (bits_remaining == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < bits_remaining && bits_remaining < P2);

        std::array<V, ARRAY_SIZE> tableVal;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller - bits_remaining);
            C cR1 = MFE_LU::getMontvalueR(mf[j]);
            tableVal[j] = MFE_LU::twoPowLimited_times_x(mf[j], index, cR1);
        }
        for (int i=0; i<bits_remaining; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   tableVal[j]);
        }
        return result;
    }
  }


#ifdef __clang__
#  pragma GCC diagnostic pop
#endif
};


}} // end namespace

#endif
