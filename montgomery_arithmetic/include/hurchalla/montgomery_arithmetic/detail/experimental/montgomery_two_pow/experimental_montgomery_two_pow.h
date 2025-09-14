// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_EXPERIMENTAL_MONTGOMERY_TWO_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_EXPERIMENTAL_MONTGOMERY_TWO_POW_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontgomeryFormExtensions.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/count_leading_zeros.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <type_traits>
#include <cstddef>
#include <array>
#include <utility>


// This file's purpose is to implement all ideas for fast two_pow, so they
// can be benchmarked to find what works best for a platform/compiler.


namespace hurchalla { namespace experimental {

// Implementation note: this is a modified version of the 2^k-ary exponentiation
// algorithm  ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring ),
// which uses optimizations knowing the base is exactly 2.
// When a particular implementation uses a table (as is typical for the 2^kary
// algorithm), we precalculate the even exponents as well as the normal odd
// exponents, in order to avoid two conditional branches that would exist in the
// main loop of the normal 2^k-ary algorithm.  This is particularly helpful for
// the array version of this function further below.
//
// We use a struct with static member functions to disallow ADL.

struct experimental_montgomery_two_pow {

  template <typename U>
  static constexpr int floor_log2(U x)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");
    // x > 0 required, but C++11 constexpr function won't allow this check
    //HPBC_CLOCKWORK_PRECONDITION2(x > 0);
    return (x <= 1) ? 0 : 1 + experimental_montgomery_two_pow::floor_log2(x >> 1);
  }


  // Calculate pow(2, n), modulo the modulus of mf, and return the result in
  // montgomeryform representation.
  template <class MF, typename U,
            bool USE_SLIDING_WINDOW_OPTIMIZATION = true,
            size_t TABLE_BITS = 5, size_t CODE_SECTION = 0,
            bool USE_SQUARING_VALUE_OPTIMIZATION = false>
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    // FYI you almost certainly want either TABLE_BITS == 0, or
    // 2 <= TABLE_BITS <= 6.  Anything above 9 is probably a very bad idea even
    // if it works since it would mean calculating 1024+ table entries!
    static_assert(0 <= TABLE_BITS && TABLE_BITS < 10, "");

    namespace hc = hurchalla;
    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using MFE = hc::detail::MontgomeryFormExtensions<MF, hc::LowlatencyTag>;
    using SV = typename MFE::SquaringValue;
    using std::size_t;

    constexpr int P = static_cast<int>(TABLE_BITS);

    // initialize the precalculation table for 2^k-ary pow algorithm
    static_assert(P >= 0, "");
    constexpr size_t TABLESIZE = 1u << P;
    static_assert(TABLESIZE >= 1, "");
    V table[TABLESIZE];

    // This simple 'for' loop would work fine, but we can init faster...
//    for (int i=1; i<TABLESIZE; ++i)
//        table[i] = mf.two_times(table[i-1]);

    // Let's do optimized initializations for the different table sizes.
    if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 1) {
        // TABLESIZE of 1 means we're not using a table- no table prep needed.

        using RU = typename MFE::RU;
        constexpr int digitsRU = hc::ut_numeric_limits<RU>::digits;
        constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));
        constexpr size_t MASK = (1u << P2) - 1u;

if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        // This is almost a copy of the main code at bottom of this function,
        // but we use convertInExtended() on the fly instead of accessing a
        // table and we use P2 instead of P.
        int shift = 0;
        if (n > MASK) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > P2);
            shift = numbits - P2;
        }

        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t index = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(index <= MASK);
        RU num = static_cast<RU>(static_cast<RU>(1) << index);
        V result = MFE::convertInExtended(mf, num);
        while (shift >= P2) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > P2 && (static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                }
            }
            shift -= P2;
            index = static_cast<size_t>(n >> shift) & MASK;
            num = static_cast<RU>(static_cast<RU>(1) << index);
            V tableVal = MFE::convertInExtended(mf, num);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(P2 > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2 - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2; ++i)
                    result = mf.square(result);
            }

            result = mf.multiply(result, tableVal);
        }
        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);

        size_t tmpmask = (1u << shift) - 1u;
        index = static_cast<size_t>(n) & tmpmask;
        num = static_cast<RU>(static_cast<RU>(1) << index);
        V tableVal = MFE::convertInExtended(mf, num);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
        V result;
        if (n <= MASK) {
            size_t loindex = static_cast<size_t>(n);
            RU num = static_cast<RU>(static_cast<RU>(1) << loindex);
            result = MFE::convertInExtended(mf, num);
            return result;
        }
        RU magicValue = MFE::getMagicValue(mf);

        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits >= (P2 + 1));

        int shift = numbits - (P2 + 1);
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(tmp <= 2u*MASK + 1u);
        // Bit P2 of tmp was the leading bit, so it should always be set.
        HPBC_CLOCKWORK_ASSERT2(((tmp >> P2) & 1u) == 1u);
        size_t loindex = tmp & MASK;
        RU num = static_cast<RU>(static_cast<RU>(1) << loindex);
        result = MFE::convertInExtended_aTimesR(mf, num, magicValue);

        while (shift >= (P2 + 1)) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while ((static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                    if (shift < (P2 + 1))
                        goto break_0_1;
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= (P2 + 1));

                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                num = static_cast<RU>(static_cast<RU>(1) << loindex);
                V val1 = MFE::convertInExtended_aTimesR(mf, num, magicValue);
                HPBC_CLOCKWORK_ASSERT2(((tmp >> P2) & 1u) == 1u);
                // since the high bit is always set, we always choose
                // val1 = convertInExtended_aTimesR()

                if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                    SV sv = MFE::getSquaringValue(mf, result);
                    static_assert((P2 + 1) > 0, "");
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1) - 1; ++i)
                        sv = MFE::squareSV(mf, sv);
                    result = MFE::squareToMontgomeryValue(mf, sv);
                } else {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                        result = mf.square(result);
                }

                result = mf.multiply(result, val1);
            }
            else {
                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                num = static_cast<RU>(static_cast<RU>(1) << loindex);
                V val1 = MFE::convertInExtended_aTimesR(mf, num, magicValue);
                V val2 = MFE::convertInExtended(mf, num);

                if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                    SV sv = MFE::getSquaringValue(mf, result);
                    static_assert((P2 + 1) > 0, "");
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1) - 1; ++i)
                        sv = MFE::squareSV(mf, sv);
                    result = MFE::squareToMontgomeryValue(mf, sv);
                } else {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                        result = mf.square(result);
                }

                size_t hibit = (tmp >> P2) & 1u;
                // val1 = (hibit == 0) ? val2 : val1;
                val1.cmov(hibit == 0, val2);
                result = mf.multiply(result, val1);
            }
        }
        if (shift == 0)
            return result;

goto break_0_1;
break_0_1:

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < (P2 + 1));

        size_t tmpmask = (1u << shift) - 1u;
        size_t index = static_cast<size_t>(n) & tmpmask;
        RU num2 = static_cast<RU>(static_cast<RU>(1) << index);
        V tableVal = MFE::convertInExtended(mf, num2);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2) {
        // This is basically a copy of code section 0, except we replace
        // calls to convertInExtended() with twoPowLimited().

        int shift = 0;
        if (n > MASK) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > P2);
            shift = numbits - P2;
        }

        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t index = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(index <= MASK);
        V result = MFE::twoPowLimited(mf, index);
        while (shift >= P2) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > P2 && (static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                }
            }
            shift -= P2;
            index = static_cast<size_t>(n >> shift) & MASK;
            V tableVal = MFE::twoPowLimited(mf, index);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(P2 > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2 - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2; ++i)
                    result = mf.square(result);
            }

            result = mf.multiply(result, tableVal);
        }
        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);

        size_t tmpmask = (1u << shift) - 1u;
        index = static_cast<size_t>(n) & tmpmask;
        V tableVal = MFE::twoPowLimited(mf, index);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 3) {
        // This is basically a copy of code section 1, except we replace calls
        // to convertInExtended() with twoPowLimited(), and we replace calls to
        // convertInExtended_aTimesR() with RTimesTwoPowLimited().
        V result;
        if (n <= MASK) {
            size_t loindex = static_cast<size_t>(n);
            result = MFE::twoPowLimited(mf, loindex);
            return result;
        }
        RU magicValue = MFE::getMagicValue(mf);

        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits >= (P2 + 1));

        int shift = numbits - (P2 + 1);
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(tmp <= 2u*MASK + 1u);
        // Bit P2 of tmp was the leading bit, so it should always be set.
        HPBC_CLOCKWORK_ASSERT2(((tmp >> P2) & 1u) == 1u);
        size_t loindex = tmp & MASK;
        result = MFE::RTimesTwoPowLimited(mf, loindex, magicValue);

        while (shift >= (P2 + 1)) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while ((static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                    if (shift < (P2 + 1))
                        goto break_0_3;
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= (P2 + 1));

                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                V val1 = MFE::RTimesTwoPowLimited(mf, loindex, magicValue);
                HPBC_CLOCKWORK_ASSERT2(((tmp >> P2) & 1u) == 1u);
                // since the high bit is always set, we always choose
                // val1 = RTimesTwoPowLimited()

                if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                    SV sv = MFE::getSquaringValue(mf, result);
                    static_assert((P2 + 1) > 0, "");
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1) - 1; ++i)
                        sv = MFE::squareSV(mf, sv);
                    result = MFE::squareToMontgomeryValue(mf, sv);
                } else {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                        result = mf.square(result);
                }

                result = mf.multiply(result, val1);
            }
            else {
                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                V val1 = MFE::RTimesTwoPowLimited(mf, loindex, magicValue);
                V val2 = MFE::twoPowLimited(mf, loindex);

                if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                    SV sv = MFE::getSquaringValue(mf, result);
                    static_assert((P2 + 1) > 0, "");
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1) - 1; ++i)
                        sv = MFE::squareSV(mf, sv);
                    result = MFE::squareToMontgomeryValue(mf, sv);
                } else {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                        result = mf.square(result);
                }

                size_t hibit = (tmp >> P2) & 1u;
                // val1 = (hibit == 0) ? val2 : val1;
                val1.cmov(hibit == 0, val2);
                result = mf.multiply(result, val1);
            }
        }
        if (shift == 0)
            return result;

goto break_0_3;
break_0_3:

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < (P2 + 1));

        size_t tmpmask = (1u << shift) - 1u;
        size_t index = static_cast<size_t>(n) & tmpmask;
        V tableVal = MFE::twoPowLimited(mf, index);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 4 && CODE_SECTION <= 9) {
        // this is a version of the scalar pow from montgomery_pow.h,
        // optimized for a base of 2.
        V result = MFE::twoPowLimited(mf, static_cast<size_t>(n) & MASK);
        V base = MFE::getMontvalueR(mf);
        n = static_cast<U>(n >> P2);
   if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 4) {
        while (n > 0) {
            if (static_cast<size_t>(n) & 1u)
                result = mf.multiply(result, base);
            base = mf.square(base);
            n = static_cast<U>(n >> 1u);
        }
   } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 5) {
        V mont_one = mf.getUnityValue();
        while (n > 0) {
            V tmp = mont_one;
            tmp.cmov((static_cast<size_t>(n) & 1u), base);
            result = mf.multiply(result, tmp);
            base = mf.square(base);
            n = static_cast<U>(n >> 1u);
        }
   } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 6) {
        V mont_one = mf.getUnityValue();
        while (true) {
            V tmp = mont_one;
            tmp.cmov((static_cast<size_t>(n) & 1u), base);
            result = mf.multiply(result, tmp);
            if (n <= 1)
                break;
            base = mf.square(base);
            n = static_cast<U>(n >> 1u);
        }
   } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 7) {
        V tmp[2];
        tmp[0] = mf.getUnityValue();
        while (n > 0) {
            tmp[1] = base;
            result = mf.multiply(result, tmp[(static_cast<size_t>(n) & 1u)]);
            base = mf.square(base);
            n = static_cast<U>(n >> 1u);
        }
   } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 8) {
        // This seems to be the 'sweet' spot for the technique I use
        // in sections 7-9; though none of these attempts seem to win
        // first place.  If we use a larger tmp table than this, it
        // requires extra multiplies overall - e.g. like section 9.
        V tmp[4];
        tmp[0] = mf.getUnityValue();
        while (n > 0) {
            V baseSqrd = mf.square(base);
            tmp[3] = mf.template multiply<LowuopsTag>(baseSqrd, base);
            tmp[1] = base;
            base = mf.square(baseSqrd);
            tmp[2] = baseSqrd;
            result = mf.template multiply<LowuopsTag>(result, tmp[(static_cast<size_t>(n) & 3u)]);
            n = static_cast<U>(n >> 2u);
        }
   } else {   // CODE_SECTION 9
        using Tag1 = LowlatencyTag;
        using Tag2 = LowuopsTag;
        V tmp[8];
        tmp[0] = mf.getUnityValue();
        while (n > 0) {
            V base2 = mf.template square<Tag1>(base);
            V base4 = mf.template square<Tag1>(base2);
            tmp[1] = base;
            tmp[2] = base2;
            tmp[3] = mf.template multiply<Tag2>(base2, base);
            tmp[4] = base4;
            tmp[5] = mf.template multiply<Tag2>(base4, base);
            base = mf.template square<Tag1>(base4);
            tmp[6] = mf.template square<Tag2>(tmp[3]);
            tmp[7] = mf.template multiply<Tag2>(tmp[3], base4);

            result = mf.template multiply<Tag2>(result, tmp[(static_cast<size_t>(n) & 7u)]);
            n = static_cast<U>(n >> 3u);
        }
   }
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 10 && CODE_SECTION <= 13) {
    
    // let's try 2kary with two tables...
    {
        // make the high table size either 2, 4, 8, or 16, depending on the CODE_SECTION
        constexpr int NUMBITS_TABLE_HIGH_SIZE = static_cast<int>(CODE_SECTION) - 9;

        static_assert(1 <= NUMBITS_TABLE_HIGH_SIZE && NUMBITS_TABLE_HIGH_SIZE <= 4, "");

        constexpr size_t TABLE_HIGH_SIZE = 1u << NUMBITS_TABLE_HIGH_SIZE;
        constexpr int NUMBITS_MASKBIG = P2 + NUMBITS_TABLE_HIGH_SIZE;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;
        std::array<C, TABLE_HIGH_SIZE> table_high;

        using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;

        static_assert(2 <= TABLE_HIGH_SIZE && TABLE_HIGH_SIZE <= 16, "");

        if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 2) {
            C cR1 = MFE_LU::getMontvalueR(mf);
            table_high[0] = cR1;                                    // R^1
            table_high[1] = mf.getCanonicalValue(mf.square(cR1));   // R^2
        }
        else if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 4) {
            C cR1 = MFE_LU::getMontvalueR(mf);
            V vR2 = mf.square(cR1);
            table_high[0] = cR1;                                         // R^1
            table_high[1] = mf.getCanonicalValue(vR2);                   // R^2
            table_high[2] = mf.getCanonicalValue(mf.multiply(vR2, cR1)); // R^3
            table_high[3] = mf.getCanonicalValue(mf.square(vR2));        // R^4
        }
        else if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 8) {
            C cR1 = MFE_LU::getMontvalueR(mf);
            V vR2 = mf.square(cR1);
            V vR3 = mf.multiply(vR2, cR1);
            V vR4 = mf.square(vR2);

            table_high[0] = cR1;                                                               // R^1
            table_high[1] = mf.getCanonicalValue(vR2);                                         // R^2
            table_high[2] = mf.getCanonicalValue(vR3);                                         // R^3
            table_high[3] = mf.getCanonicalValue(vR4);                                         // R^4

            table_high[4] = mf.getCanonicalValue(mf.template multiply<LowuopsTag>(vR3, vR2));  // R^5
            table_high[5] = mf.getCanonicalValue(mf.template square<LowuopsTag>(vR3));         // R^6
            table_high[6] = mf.getCanonicalValue(mf.multiply(vR4, vR3));                       // R^7
            table_high[7] = mf.getCanonicalValue(mf.square(vR4));                              // R^8
        }
        else if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 16) {
            C cR1 = MFE_LU::getMontvalueR(mf);
            V vR2 = mf.square(cR1);
            V vR3 = mf.multiply(vR2, cR1);
            V vR4 = mf.square(vR2);
            table_high[0] = cR1;                                                               // R^1
            table_high[1] = mf.getCanonicalValue(vR2);                                         // R^2
            table_high[2] = mf.getCanonicalValue(vR3);                                         // R^3
            table_high[3] = mf.getCanonicalValue(vR4);                                         // R^4

            V vR5 = mf.multiply(vR3, vR2);
            table_high[4] = mf.getCanonicalValue(vR5);                                         // R^5
            table_high[5] = mf.getCanonicalValue(mf.square(vR3));                              // R^6
            table_high[6] = mf.getCanonicalValue(mf.template multiply<LowuopsTag>(vR4, vR3));  // R^7
            table_high[7] = mf.getCanonicalValue(mf.template square<LowuopsTag>(vR4));         // R^8

            table_high[8] = mf.getCanonicalValue(mf.template multiply<LowuopsTag>(vR5, vR4));                       // R^9
            table_high[9] = mf.getCanonicalValue(mf.template square<LowuopsTag>(vR5));                              // R^10
            table_high[10] = mf.getCanonicalValue(mf.template multiply<LowuopsTag>(table_high[5], vR5));            // R^11
            table_high[11] = mf.getCanonicalValue(mf.template square<LowuopsTag>(table_high[5]));                   // R^12
            table_high[12] = mf.getCanonicalValue(mf.template multiply<LowuopsTag>(table_high[6], table_high[5]));  // R^13
            table_high[13] = mf.getCanonicalValue(mf.template square<LowuopsTag>(table_high[6]));                   // R^14
            table_high[14] = mf.getCanonicalValue(mf.multiply(table_high[7], table_high[6]));                       // R^15
            table_high[15] = mf.getCanonicalValue(mf.square(table_high[7]));                                        // R^16
        }

        int shift = 0;
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
        }

        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        size_t loindex = tmp & MASK;
        size_t hiindex = tmp >> P2;
        HPBC_CLOCKWORK_ASSERT2(hiindex < TABLE_HIGH_SIZE);
        V result = MFE::twoPowLimited_times_x(mf, loindex, table_high[hiindex]);

        while (shift >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > NUMBITS_MASKBIG && (static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                }
            }
            HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

            shift -= NUMBITS_MASKBIG;
            tmp = static_cast<size_t>(n >> shift);
            loindex = tmp & MASK;
            hiindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, table_high[hiindex]);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(NUMBITS_MASKBIG > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);
            }

            result = mf.multiply(result, val1);
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        size_t index = static_cast<size_t>(n) & tmpmask;
        loindex = index & MASK;
        hiindex = (index >> P2) & (TABLE_HIGH_SIZE - 1);
        V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, table_high[hiindex]);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val1);
        return result;
    }
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 14 && CODE_SECTION <= 16) {
    // why stop at two tables?  Now 2kary with three tables!

    using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;

    static_assert(CODE_SECTION >= 14 && CODE_SECTION <= 16, "");
    {
        // make the high tables' size either 2, 4, 8, depending on the CODE_SECTION
        constexpr int NUMBITS_TABLE_HIGH_SIZE = static_cast<int>(CODE_SECTION) - 13;

        static_assert(1 <= NUMBITS_TABLE_HIGH_SIZE && NUMBITS_TABLE_HIGH_SIZE <= 3, "");

        constexpr size_t TABLE_HIGH_SIZE = 1u << NUMBITS_TABLE_HIGH_SIZE;
        constexpr int NUMBITS_MASKBIG = P2 + NUMBITS_TABLE_HIGH_SIZE + NUMBITS_TABLE_HIGH_SIZE;
        constexpr int P3 = P2 + NUMBITS_TABLE_HIGH_SIZE;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;
        std::array<C, TABLE_HIGH_SIZE> table_mid;
        std::array<V, TABLE_HIGH_SIZE> table_high;

        static_assert(2 <= TABLE_HIGH_SIZE && TABLE_HIGH_SIZE <= 8, "");
        if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE == 2) {
            C r1 = MFE::getMontvalueR(mf);
            V r2 = mf.square(r1);
            table_mid[0] = r1;                         // R^1
            table_mid[1] = mf.getCanonicalValue(r2);   // R^2

            table_high[0] = mf.getUnityValue();   // R^0
            table_high[1] = r2;                   // R^2
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

            V r8 = mf.square(r4);
            table_high[0] = mf.getUnityValue();   // R^0
            table_high[1] = r4;                   // R^4
            table_high[2] = r8;                   // R^8
            table_high[3] = mf.multiply(r8, r4);  // R^12
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

            V r24 = mf.multiply(r16, r8);
            V r32 = mf.square(r16);

            table_high[0] = mf.getUnityValue();                          // R^0
            table_high[1] = r8;                                          // R^8
            table_high[2] = r16;                                         // R^16
            table_high[3] = r24;                                         // R^24
            table_high[4] = r32;                                         // R^32
            table_high[5] = mf.template multiply<LowuopsTag>(r24, r16);  // R^40
            table_high[6] = mf.square(r24);                              // R^48
            table_high[7] = mf.multiply(r32, r24);                       // R^56
        }

        int shift = 0;
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
        }

        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        size_t loindex = tmp & MASK;
        size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
        V ttx = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

        size_t hiindex = tmp >> P3;
        HPBC_CLOCKWORK_ASSERT2(hiindex < TABLE_HIGH_SIZE);
        V result = mf.multiply(ttx, table_high[hiindex]);

        while (shift >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > NUMBITS_MASKBIG && (static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                }
            }
            HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

            shift -= NUMBITS_MASKBIG;
            tmp = static_cast<size_t>(n >> shift);
            loindex = tmp & MASK;
            midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            hiindex = (tmp >> P3) & (TABLE_HIGH_SIZE - 1);

            ttx = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                static_assert(NUMBITS_MASKBIG >= 3, "");
                SV sv = MFE::getSquaringValue(mf, result);
                sv = MFE::squareSV(mf, sv);
                sv = MFE::squareSV(mf, sv);

                V val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=2; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);

                result = mf.multiply(result, val1);
            }
            else {
                static_assert(NUMBITS_MASKBIG >= 2, "");
                result = mf.square(result);
                result = mf.square(result);

                V val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=2; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);

                result = mf.multiply(result, val1);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        size_t index = static_cast<size_t>(n) & tmpmask;
        loindex = index & MASK;
        midindex = (index >> P2) & (TABLE_HIGH_SIZE - 1);
        ttx = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

        hiindex = (index >> P3) & (TABLE_HIGH_SIZE - 1);
        V val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val1);
        return result;
    }
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 17) {
        // This is an updated version of CODE_SECTION 3 to use
        // the new and more general MFE functions.

        constexpr int NUMBITS_EXTRA = 1;

        constexpr int NUMBITS_MASKBIG = P2 + NUMBITS_EXTRA;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

        using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;
        C cR1 = MFE::getMontvalueR(mf);
        C cR2 = mf.getCanonicalValue(mf.template square<LowlatencyTag>(cR1));

        V result;
        if (n <= MASKBIG) {
            size_t loindex = static_cast<size_t>(n) & MASK;
            C cHigh = cR1;
            cHigh.cmov(((static_cast<size_t>(n) >> (NUMBITS_MASKBIG - 1)) & 1u), cR2);
            result = MFE::twoPowLimited_times_x(mf, loindex, cHigh);
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);

        int shift = numbits - NUMBITS_MASKBIG;
        HPBC_CLOCKWORK_ASSERT2(shift > 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        // we know the leading bit of n is by definition set.
        HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 1u);
#if 0 
        // due to above assert, we don't need this section
        C cHigh = cR1;
        cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR2);
#else
        C cHigh = cR2;
#endif
        size_t loindex = tmp & MASK;
        result = MFE::twoPowLimited_times_x(mf, loindex, cHigh);

        while (shift >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while ((static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                    if (shift < NUMBITS_MASKBIG)
                        goto break_0_17;
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

                shift -= NUMBITS_MASKBIG;
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 1u);
                // since the high bit is always set, we always choose
                cHigh = cR2;
            }
            else {
                shift -= NUMBITS_MASKBIG;
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                cHigh = cR1;
                cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR2);
            }

            V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, cHigh);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(NUMBITS_MASKBIG > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);
            }

            result = mf.multiply(result, val1);
        }
        if (shift == 0)
            return result;

goto break_0_17;
break_0_17:

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        tmp = static_cast<size_t>(n) & tmpmask;
#if 0
        // this would be the general solution, but our CODE_SECTION here is a
        // special case - we can use the #else instead.
        loindex = tmp & MASK;
        cHigh = cR1;
        cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR2);
#else
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
        loindex = tmp;
        HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 0u);
        cHigh = cR1;
#endif
        V val1 = MFE::twoPowLimited_times_x(mf, loindex, cHigh);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val1);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 18) {
        // This extends CODE_SECTION 17 to use 2 extra bits instead of 1

        constexpr int NUMBITS_EXTRA = 2;

        constexpr int NUMBITS_MASKBIG = P2 + NUMBITS_EXTRA;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;
        using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;

        C cR1 = MFE::getMontvalueR(mf);
        V vR2 = mf.square(cR1);
        V vR3 = mf.multiply(vR2, cR1);
        V vR4 = mf.square(vR2);

        C cR2 = mf.getCanonicalValue(vR2);
        C cR3 = mf.getCanonicalValue(vR3);
        C cR4 = mf.getCanonicalValue(vR4);

        V result;
        if (n <= MASKBIG) {
            size_t tmp = static_cast<size_t>(n);
            size_t loindex = tmp & MASK;
            C cHighx0 = cR1;
            cHighx0.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR3);
            C cHighx1 = cR2;
            cHighx1.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR4);
            C cHigh = cHighx0;
            cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 2)) & 1u), cHighx1);
            result = MFE::twoPowLimited_times_x(mf, loindex, cHigh);
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);

        int shift = numbits - NUMBITS_MASKBIG;
        HPBC_CLOCKWORK_ASSERT2(shift > 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);

        // we know the leading bit of n is by definition set.
        HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 1u);
        C cHigh = cR3;
        cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 2)) & 1u), cR4);

        size_t loindex = tmp & MASK;
        result = MFE::twoPowLimited_times_x(mf, loindex, cHigh);

        while (shift >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while ((static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                    if (shift < NUMBITS_MASKBIG)
                        goto break_0_18;
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

                shift -= NUMBITS_MASKBIG;
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 1u);
                cHigh = cR3;
                cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 2)) & 1u), cR4);
            }
            else {
                shift -= NUMBITS_MASKBIG;
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                C cHighx0 = cR1;
                cHighx0.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR3);
                C cHighx1 = cR2;
                cHighx1.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR4);
                cHigh = cHighx0;
                cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 2)) & 1u), cHighx1);
            }

            V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, cHigh);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(NUMBITS_MASKBIG > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);
            }

            result = mf.multiply(result, val1);
        }
        if (shift == 0)
            return result;

goto break_0_18;
break_0_18:

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        tmp = static_cast<size_t>(n) & tmpmask;
        loindex = tmp & MASK;

        HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 0u);
        cHigh = cR1;
        cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 2)) & 1u), cR2);

        V val1 = MFE::twoPowLimited_times_x(mf, loindex, cHigh);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val1);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 19 && CODE_SECTION <= 21) {
    // If three tables is good, then four is obviously better(?)
    //
    // this is an expanded version of CODE_SECTIONS 14 to 16, using 4 tables instead of a puny 3

    static_assert(CODE_SECTION >= 19 && CODE_SECTION <= 21, "");

    using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;

    {
        // make the high tables' size either 2, 4, 8, depending on the CODE_SECTION
        constexpr int NUMBITS_TABLE_HIGH_SIZE = static_cast<int>(CODE_SECTION) - 18;

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
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
        }

        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(n >> shift);
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


        while (shift >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > NUMBITS_MASKBIG && (static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                }
            }
            HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

            shift -= NUMBITS_MASKBIG;
            tmp = static_cast<size_t>(n >> shift);
            loindex = tmp & MASK;
            midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            hiindex = (tmp >> P3) & (TABLE_HIGH_SIZE - 1);
            ultindex = (tmp >> P4) & (TABLE_HIGH_SIZE - 1);

            V ttx = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            V val1;
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                static_assert(NUMBITS_MASKBIG >= 5, "");
                SV sv = MFE::getSquaringValue(mf, result);
                sv = MFE::squareSV(mf, sv);
                sv = MFE::squareSV(mf, sv);

                val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

                sv = MFE::squareSV(mf, sv);
                sv = MFE::squareSV(mf, sv);

                val1 = mf.template multiply<LowuopsTag>(val1, table_ultra[ultindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=4; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            }
            else {
                static_assert(NUMBITS_MASKBIG >= 4, "");
                result = mf.square(result);
                result = mf.square(result);

                val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);

                result = mf.square(result);
                result = mf.square(result);

                val1 = mf.template multiply<LowuopsTag>(val1, table_ultra[ultindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=4; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);
            }

            result = mf.multiply(result, val1);
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        size_t index = static_cast<size_t>(n) & tmpmask;
        loindex = index & MASK;
        midindex = (index >> P2) & (TABLE_HIGH_SIZE - 1);
        hiindex = (index >> P3) & (TABLE_HIGH_SIZE - 1);
        ultindex = (index >> P4) & (TABLE_HIGH_SIZE - 1);
        V ttx = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

        result = mf.square(result);

        V val1 = mf.template multiply<LowuopsTag>(ttx, table_high[hiindex]);
        val1 = mf.template multiply<LowuopsTag>(val1, table_ultra[ultindex]);

        for (int i=1; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val1);
        return result;
    }
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 22 && CODE_SECTION <= 26) {
// super duper experimental.  lots of extra tables, all size 4.

    using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;
    {
        constexpr int NUMBITS_TABLE_HIGH_SIZE = 2;
        constexpr int NUM_EXTRA_TABLES = 2 * (CODE_SECTION - 21);

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
            size_t tmp = static_cast<size_t>(n >> shift);
            HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
            size_t loindex = tmp & MASK;
            size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            result = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            V next = r4;              // R^4
            HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i < NUM_EXTRA_TABLES; ++i) {
                tables_extra[i][0] = mf.getUnityValue();   // R^0
                tables_extra[i][1] = next;
                V nextSq = mf.square(next);
                V nexttmp = mf.square(nextSq);
                tables_extra[i][2] = nextSq;
                tables_extra[i][3] = mf.template multiply<LowuopsTag>(nextSq, next);
                next = nexttmp;

                int P_extra = P3 + i * NUMBITS_TABLE_HIGH_SIZE;
                size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                result = mf.template multiply<LowuopsTag>(tables_extra[i][index_extra], result);
            }
        }


        while (shift >= NUMBITS_MASKBIG) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > NUMBITS_MASKBIG && (static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                }
            }
            HPBC_CLOCKWORK_ASSERT2(shift >= NUMBITS_MASKBIG);

            shift -= NUMBITS_MASKBIG;
            size_t tmp = static_cast<size_t>(n >> shift);
            size_t loindex = tmp & MASK;
            size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            V val1 = MFE_LU::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P3-1; ++i)
                    sv = MFE::squareSV(mf, sv);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i < NUM_EXTRA_TABLES; ++i) {
                    int P_extra = P3 + i * NUMBITS_TABLE_HIGH_SIZE;
                    size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                    val1 = mf.template multiply<LowuopsTag>(val1, tables_extra[i][index_extra]);

                    static_assert(NUMBITS_TABLE_HIGH_SIZE == 2, "");
                    sv = MFE::squareSV(mf, sv);
                    sv = MFE::squareSV(mf, sv);
                }
                result = MFE::squareToMontgomeryValue(mf, sv);
            }
            else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P3; ++i)
                    result = mf.square(result);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i < NUM_EXTRA_TABLES; ++i) {
                    int P_extra = P3 + i * NUMBITS_TABLE_HIGH_SIZE;
                    size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
                    val1 = mf.template multiply<LowuopsTag>(val1, tables_extra[i][index_extra]);

                    static_assert(NUMBITS_TABLE_HIGH_SIZE == 2, "");
                    result = mf.square(result);
                    result = mf.square(result);
                }
            }

            result = mf.multiply(result, val1);
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        size_t tmp = static_cast<size_t>(n) & tmpmask;
        size_t loindex = tmp & MASK;
        size_t midindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
        V val1 = MFE::twoPowLimited_times_x(mf, loindex, table_mid[midindex]);

        result = mf.square(result);

        // could use:
        //for (int i=0; i < NUM_EXTRA_TABLES && (2*i + P3 < shift); ++i)

        HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i < NUM_EXTRA_TABLES; ++i) {
            int P_extra = P3 + i * NUMBITS_TABLE_HIGH_SIZE;
            size_t index_extra = (tmp >> P_extra) & (TABLE_HIGH_SIZE - 1);
            val1 = mf.multiply(val1, tables_extra[i][index_extra]);
        }

        for (int i=1; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val1);
        return result;
    }
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 27) {
// This is a quite elegant optimization of CODE_SECTION 2, with very low uops.
// The real speed of this should come when adapted for the array two_pow, due to
// the extremely low uops.

        int shift = 0;
        if (n > MASK) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > P2);
            shift = numbits - P2;
        }

        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t index = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(index <= MASK);
        V result = MFE::twoPowLimited(mf, index);

        while (shift >= P2) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > P2 && (static_cast<size_t>(n>>(shift-1)) & 1u) == 0) {
                    result = mf.square(result);
                    --shift;
                }
            }

            // Multiplying by 2 here becomes a multiply by R after P2 squarings.
            // At the end of this loop iteration, the extra factor R will be
            // removed from result by the twoPowLimited_times_x() call (which
            // requires 'x' to have an extra factor of R).
            result = mf.two_times(result);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(P2 > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2 - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2; ++i)
                    result = mf.square(result);
            }

            shift -= P2;
            index = static_cast<size_t>(n >> shift) & MASK;
            C tmp = mf.getCanonicalValue(result);
            result = MFE::twoPowLimited_times_x(mf, index, tmp);
        }
        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);

        size_t tmpmask = (1u << shift) - 1u;
        index = static_cast<size_t>(n) & tmpmask;
        V tableVal = MFE::twoPowLimited(mf, index);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 28) {
// This is a further optimized version of CODE_SECTION 27 for low uops.

        int shift = 0;
        if (n > MASK) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > P2);
            shift = numbits - P2;
        }

        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t index = static_cast<size_t>(n >> shift);
        HPBC_CLOCKWORK_ASSERT2(index <= MASK);
        C cR1 = MFE::getMontvalueR(mf);
        V result = MFE::twoPowLimited_times_x_v2(mf, index + 1, cR1);

        while (shift >= P2) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(P2 > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2 - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2; ++i)
                    result = mf.square(result);
            }

            shift -= P2;
            index = static_cast<size_t>(n >> shift) & MASK;
            C tmp = mf.getCanonicalValue(result);
            result = MFE::twoPowLimited_times_x_v2(mf, index + 1, tmp);
        }
        result = mf.divideBySmallPowerOf2(mf.getCanonicalValue(result), 1);

        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);

        size_t tmpmask = (1u << shift) - 1u;
        index = static_cast<size_t>(n) & tmpmask;
        V tableVal = MFE::twoPowLimited(mf, index);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 29) {
// This is a further optimized version of CODE_SECTION 28 for low uops.

        int shift = 0;
        if (n > MASK) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > P2);
            shift = numbits - P2;
        }

        C cR1 = MFE::getMontvalueR(mf);
        C cresult = cR1;

        while (shift >= P2) {
            size_t index = static_cast<size_t>(n >> shift) & MASK;
            V result = MFE::twoPowLimited_times_x_v2(mf, index + 1, cresult);

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                static_assert(P2 > 0, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2 - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            } else {
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2; ++i)
                    result = mf.square(result);
            }
            cresult = mf.getCanonicalValue(result);

            shift -= P2;
        }
        size_t index = static_cast<size_t>(n >> shift) & MASK;
        V result = MFE::twoPowLimited_times_x(mf, index, cresult);

        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);

        size_t tmpmask = (1u << shift) - 1u;
        index = static_cast<size_t>(n) & tmpmask;
        V tableVal = MFE::twoPowLimited_times_x(mf, index, cR1);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
} else {
    static_assert(CODE_SECTION == 30, "");
    // for comparison purposes, this is the current MontgomeryForm pow.
    // the static cast may lose bits, so this might not be an exact benchmark
    C mont_two = mf.two_times(mf.getUnityValue());
    return mf.pow(mont_two, static_cast<typename MF::IntegerType>(n));
}
    }
    else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 2) {
        table[0] = mf.getUnityValue();   // montgomery one

// The different code sections should be functionally equivalent.  You can test
// to see which is fastest.
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        table[1] = mf.two_times(table[0]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1 || CODE_SECTION == 2) {
// For this particular case of TABLE_BITS == 1 (TABLESIZE == 2), we can
// use a version of 2^k-ary that is heavily optimized for the 1 bit table:
        C mont_one = mf.getUnityValue();
        C mont_two = mf.two_times(mont_one);
        V result;
        if (n <= 1) {
            result = (n == 0) ? mont_one : mont_two;
            return result;
        }
        HPBC_CLOCKWORK_ASSERT2(n > 1);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > 1);
        int shift = numbits - 1;
        HPBC_CLOCKWORK_ASSERT2((n >> shift) == 1);   // because shift == numbits - 1
        C cresult = mont_two;

        HPBC_CLOCKWORK_ASSERT2(shift >= 1);
        // since cresult == 2, this two_times() is equivalent to squaring
        cresult = mf.two_times(cresult);
        --shift;
  if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {   // branch based code
        if (static_cast<size_t>(n >> shift) & 1)
            cresult = mf.two_times(cresult);
  } else {   // cmov based code
        C ctmp = mf.two_times(cresult);
        cresult.cmov(static_cast<size_t>(n >> shift) & 1, ctmp);
  }
        result = cresult;

        while (shift >= 1) {
            result = mf.square(result);
            --shift;
  if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {   // branch based code
            if (static_cast<size_t>(n >> shift) & 1)
                result = mf.two_times(result);
  } else {   // cmov based code
            V vtmp = mf.two_times(result);
            result.cmov(static_cast<size_t>(n >> shift) & 1, vtmp);
  }
        }
        return result;
} else {
        static_assert(CODE_SECTION == 3 || CODE_SECTION == 4, "");
        // This is a better optimized version of (just above) CODE_SECTIONs 1, 2
        int shift;
        V result;
        {
            // This portion is basically a copy/paste of the setup code from
            // TABLE_BITS 0's CODE_SECTION 17.  It's likely the fastest setup we
            // can use for this code section.
            // ------
            using RU = typename MFE::RU;
            constexpr int digitsRU = hc::ut_numeric_limits<RU>::digits;
            constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));
            constexpr size_t MASK = (1u << P2) - 1u;

            constexpr int NUMBITS_MASKBIG = P2 + 1;
            constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

            C cR1 = MFE::getMontvalueR(mf);
            C cR2 = mf.getCanonicalValue(mf.template square<LowlatencyTag>(cR1));

            if (n <= MASKBIG) {
                size_t loindex = static_cast<size_t>(n) & MASK;
                C cHigh = cR1;
                cHigh.cmov(((static_cast<size_t>(n) >> (NUMBITS_MASKBIG - 1)) & 1u), cR2);
                result = MFE::twoPowLimited_times_x(mf, loindex, cHigh);
                return result;
            }

            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);

            shift = numbits - NUMBITS_MASKBIG;
            HPBC_CLOCKWORK_ASSERT2(shift > 0);
            size_t tmp = static_cast<size_t>(n >> shift);
            HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
            // we know the leading bit of n is by definition set.
            HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 1u);
            size_t loindex = tmp & MASK;
            result = MFE::twoPowLimited_times_x(mf, loindex, cR2);
        }

        while (shift >= 1) {
            result = mf.square(result);
            --shift;
  if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 3) {   // branch based code
            if (static_cast<size_t>(n >> shift) & 1)
                result = mf.two_times(result);
  } else {   // cmov based code
            V vtmp = mf.two_times(result);
            result.cmov(static_cast<size_t>(n >> shift) & 1, vtmp);
  }
        }
        return result;
}
    } else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 4) {
        C m1 = mf.getUnityValue();   // montgomery one
        C m2 = mf.two_times(m1);
        C m4 = mf.two_times(m2);
        C m8 = mf.two_times(m4);
        table[0] = m1;
        table[1] = m2;
        table[2] = m4;
        table[3] = m8;

if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
} else {
// try the same optimization as TABLE_BITS == 1 (TABLESIZE == 2), but
// using the 4 entry table that we have here.

        // this section is a copy/paste of the main code at bottom of this
        // function.  Note P == 2 here.
        //
        V result;
        constexpr size_t MASK = TABLESIZE - 1;
        if (n <= MASK) {
            result = table[static_cast<size_t>(n)];
            return result;
        }
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P);
        int shift = numbits - P;
        U tmp = n >> shift;
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
        size_t index = static_cast<size_t>(tmp);
        result = table[index];

        // this section is not copy/paste
        //
        while (shift >= 1) {
            result = mf.square(result);
            --shift;
  if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
            if (static_cast<size_t>(n >> shift) & 1)
                result = mf.two_times(result);
  } else {
            V vtmp = mf.two_times(result);
            result.cmov(static_cast<size_t>(n >> shift) & 1, vtmp);
  }
        }
        return result;
}
    } else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 8) {
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        // for MontyHalfRange and MontyFullRangeMasked, this would be slightly
        // faster with a table of C (CanonicalValue) rather than table's type
        // which is V (MontgomeryValue).  Since I don't expect this code section
        // to ever be competitive with the fastest methods, I'm not expecting to
        // micro-optimize it more with a table of type C.
        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);
        table[7] = mf.two_times(table[6]);
} else {
// try the same optimization as TABLE_BITS == 1 (TABLESIZE == 2), but
// using the 8 entry table that we have here.

        C ctable[TABLESIZE];
        ctable[0] = mf.getUnityValue();   // montgomery one
        ctable[1] = mf.two_times(ctable[0]);
        ctable[2] = mf.two_times(ctable[1]);
        ctable[3] = mf.two_times(ctable[2]);
        ctable[4] = mf.two_times(ctable[3]);
        ctable[5] = mf.two_times(ctable[4]);
        ctable[6] = mf.two_times(ctable[5]);
        ctable[7] = mf.two_times(ctable[6]);
        // this section is a copy/paste of the main code at bottom of this
        // function.  Note P == 3 here.
        //
        V result;
        constexpr size_t MASK = TABLESIZE - 1;
        if (n <= MASK) {
            result = ctable[static_cast<size_t>(n)];
            return result;
        }
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P);
        int shift = numbits - P;
        U tmp = n >> shift;
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
        size_t index = static_cast<size_t>(tmp);
        result = ctable[index];

        // this section is not copy/paste
        //
        while (shift >= 1) {
            result = mf.square(result);
            --shift;
  if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
            if (static_cast<size_t>(n >> shift) & 1)
                result = mf.two_times(result);
  } else {
            V vtmp = mf.two_times(result);
            result.cmov(static_cast<size_t>(n >> shift) & 1, vtmp);
  }
        }
        return result;
}
    } else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 16) {
        // for MontyHalfRange and MontyFullRangeMasked, the two_times calls
        // would be slightly faster if table were type C (and using
        // getCanonicalValue as needed), but I don't expect these code sections
        // to ever be fast enough to be worth optimizing with table type C.

if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        table[0] = mf.getUnityValue();   // montgomery one
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<15; ++i)
            table[i+1] = mf.two_times(table[i]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
        static_assert(11 < ut_numeric_limits<typename MFE::RU>::digits, "");
        table[11] = MFE::twoPowLimited(mf, 11);

        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);

        table[7] = mf.two_times(table[6]);
        table[12] = mf.two_times(table[11]);

        table[8]  = mf.two_times(table[7]);
        table[13] = mf.two_times(table[12]);

        table[9]  = mf.two_times(table[8]);
        table[14] = mf.two_times(table[13]);

        table[10] = mf.two_times(table[9]);
        table[15] = mf.two_times(table[14]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2) {
        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);

        table[12] = mf.square(table[6]);
        table[7] = mf.two_times(table[6]);
        table[14] = mf.square(table[7]);
        table[8]  = mf.two_times(table[7]);
        table[9]  = mf.two_times(table[8]);
        table[10] = mf.two_times(table[9]);

        table[11] = mf.two_times(table[10]);
        table[13] = mf.two_times(table[12]);

        table[15] = mf.two_times(table[14]);
} else {
        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);
        table[7] = mf.two_times(table[6]);
        table[14] = mf.square(table[7]);

        table[8]  = mf.two_times(table[7]);
        table[9]  = mf.two_times(table[8]);
        table[10] = mf.two_times(table[9]);
        table[11] = mf.two_times(table[10]);
        table[12] = mf.two_times(table[11]);

        table[13] = mf.two_times(table[12]);
        table[15] = mf.two_times(table[14]);
}
    }
    else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 32) {

        // for MontyHalfRange and MontyFullRangeMasked, the two_times calls
        // would be slightly faster if table were type C (and using
        // getCanonicalValue as needed), but I don't expect these code sections
        // to ever be fast enough to be worth optimizing with table type C.
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        table[0] = mf.getUnityValue();   // montgomery one
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<31; ++i)
            table[i+1] = mf.two_times(table[i]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
//        note: if using multiply, it would have to be moved lower
//        table[19] = mf.multiply(table[9]), table[10]);

        static_assert(19 < ut_numeric_limits<typename MFE::RU>::digits, "");
        table[19] = MFE::twoPowLimited(mf, 19);

        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);

        table[7] = mf.two_times(table[6]);
        table[20] = mf.two_times(table[19]);

        table[8] = mf.two_times(table[7]);
        table[21] = mf.two_times(table[20]);

        table[9]  = mf.two_times(table[8]);
        table[22] = mf.two_times(table[21]);

        table[10] = mf.two_times(table[9]);
        table[23] = mf.two_times(table[22]);

        table[11] = mf.two_times(table[10]);
        table[24] = mf.two_times(table[23]);

        table[12] = mf.two_times(table[11]);
        table[25] = mf.two_times(table[24]);

        table[13] = mf.two_times(table[12]);
        table[26] = mf.two_times(table[25]);

        table[14] = mf.two_times(table[13]);
        table[27] = mf.two_times(table[26]);

        table[15] = mf.two_times(table[14]);
        table[28] = mf.two_times(table[27]);

        table[16] = mf.two_times(table[15]);
        table[29] = mf.two_times(table[28]);

        table[17] = mf.two_times(table[16]);
        table[30] = mf.two_times(table[29]);

        table[18] = mf.two_times(table[17]);
        table[31] = mf.two_times(table[30]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2) {
//        note: if using squares, they would have to be moved lower
//        table[16] = mf.square(table[8]);
//        table[24] = mf.square(table[12]);

        static_assert(16 < ut_numeric_limits<typename MFE::RU>::digits, "");
        table[16] = MFE::twoPowLimited(mf, 16);
        static_assert(24 < ut_numeric_limits<typename MFE::RU>::digits, "");
        table[24] = MFE::twoPowLimited(mf, 24);

        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);
        table[7] = mf.two_times(table[6]);
        table[8] = mf.two_times(table[7]);

        table[9]  = mf.two_times(table[8]);
        table[17] = mf.two_times(table[16]);
        table[25] = mf.two_times(table[24]);

        table[10] = mf.two_times(table[9]);
        table[18] = mf.two_times(table[17]);
        table[26] = mf.two_times(table[25]);

        table[11] = mf.two_times(table[10]);
        table[19] = mf.two_times(table[18]);
        table[27] = mf.two_times(table[26]);

        table[12] = mf.two_times(table[11]);
        table[20] = mf.two_times(table[19]);
        table[28] = mf.two_times(table[27]);

        table[13] = mf.two_times(table[12]);
        table[21] = mf.two_times(table[20]);
        table[29] = mf.two_times(table[28]);

        table[14] = mf.two_times(table[13]);
        table[22] = mf.two_times(table[21]);
        table[30] = mf.two_times(table[29]);

        table[15] = mf.two_times(table[14]);
        table[23] = mf.two_times(table[22]);
        table[31] = mf.two_times(table[30]);
} else {
        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);
        table[7] = mf.two_times(table[6]);
        table[8] = mf.two_times(table[7]);

    // Note: for some of the clauses below, in principle we could call
    // convertIn() - note doing so would tie us to using a MontgomeryForm
    // that computes RSquaredModN in its constructor.  [At this point it
    // seems we'll probably want to compute RSquaredModN in general though,
    // to take advantage of a big speed up from twoPowLimited().  So it
    // probably wouldn't be any problem to use convertIn(), or better yet
    // to use twoPowLimited().]
  if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 3) {
        table[9]  = mf.two_times(table[8]);
        table[10] = mf.two_times(table[9]);
        table[11] = mf.two_times(table[10]);
        table[12] = mf.two_times(table[11]);
        table[13] = mf.two_times(table[12]);

        table[26] = mf.square(table[13]);
//        table[26] = mf.convertIn(4 * 256 * 256 * 256);
        table[14] = mf.two_times(table[13]);
        table[15] = mf.two_times(table[14]);
        table[16] = mf.two_times(table[15]);
        table[17] = mf.two_times(table[16]);
        table[18] = mf.two_times(table[17]);
        table[19] = mf.two_times(table[18]);
        table[20] = mf.two_times(table[19]);

        table[21] = mf.two_times(table[20]);
        table[27] = mf.two_times(table[26]);

        table[22] = mf.two_times(table[21]);
        table[28] = mf.two_times(table[27]);

        table[23] = mf.two_times(table[22]);
        table[29] = mf.two_times(table[28]);

        table[24] = mf.two_times(table[23]);
        table[30] = mf.two_times(table[29]);

        table[25] = mf.two_times(table[24]);
        table[31] = mf.two_times(table[30]);
  } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 4) {
        table[9]  = mf.two_times(table[8]);

        table[18] = mf.square(table[9]);
        table[10] = mf.two_times(table[9]);
        table[11] = mf.two_times(table[10]);
        table[12] = mf.two_times(table[11]);
        table[13] = mf.two_times(table[12]);

        table[27] = mf.multiply(table[18], table[9]);
        table[14] = mf.two_times(table[13]);
        table[19] = mf.two_times(table[18]);
        table[15] = mf.two_times(table[14]);
        table[20] = mf.two_times(table[19]);
        table[16] = mf.two_times(table[15]);
        table[21] = mf.two_times(table[20]);
        table[17] = mf.two_times(table[16]);
        table[22] = mf.two_times(table[21]);

        table[23] = mf.two_times(table[22]);
        table[28] = mf.two_times(table[27]);
        table[24] = mf.two_times(table[23]);
        table[29] = mf.two_times(table[28]);
        table[25] = mf.two_times(table[24]);
        table[30] = mf.two_times(table[29]);
        table[26] = mf.two_times(table[25]);
        table[31] = mf.two_times(table[30]);
  } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 5) {
        table[16] = mf.square(table[8]);
        table[9]  = mf.two_times(table[8]);
        table[10] = mf.two_times(table[9]);
        table[11] = mf.two_times(table[10]);
        table[12] = mf.two_times(table[11]);

        table[24] = mf.square(table[12]);
        table[28] = mf.multiply(table[12], table[16]);

        table[13] = mf.two_times(table[12]);
        table[17] = mf.two_times(table[16]);
        table[14] = mf.two_times(table[13]);
        table[18] = mf.two_times(table[17]);
        table[15] = mf.two_times(table[14]);
        table[19] = mf.two_times(table[18]);

        table[20] = mf.two_times(table[19]);

        table[21] = mf.two_times(table[20]);
        table[25] = mf.two_times(table[24]);
        table[29] = mf.two_times(table[28]);

        table[22] = mf.two_times(table[21]);
        table[26] = mf.two_times(table[25]);
        table[30] = mf.two_times(table[29]);

        table[23] = mf.two_times(table[22]);
        table[27] = mf.two_times(table[26]);
        table[31] = mf.two_times(table[30]);
  } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 6) {
        table[9]  = mf.two_times(table[8]);
        table[10] = mf.two_times(table[9]);
        table[11] = mf.two_times(table[10]);

        table[22] = mf.square(table[11]);

        table[12] = mf.two_times(table[11]);
        table[13] = mf.two_times(table[12]);
        table[14] = mf.two_times(table[13]);
        table[15] = mf.two_times(table[14]);

        table[29] = mf.multiply(table[15], table[14]);

        table[16] = mf.two_times(table[15]);
        table[17] = mf.two_times(table[16]);

        table[18] = mf.two_times(table[17]);
        table[23] = mf.two_times(table[22]);
        table[19] = mf.two_times(table[18]);
        table[24] = mf.two_times(table[23]);
        table[20] = mf.two_times(table[19]);
        table[25] = mf.two_times(table[24]);
        table[21] = mf.two_times(table[20]);
        table[26] = mf.two_times(table[25]);

        table[27] = mf.two_times(table[26]);
        table[30] = mf.two_times(table[29]);
        table[28] = mf.two_times(table[27]);
        table[31] = mf.two_times(table[30]);
  } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 7) {
        table[9]  = mf.two_times(table[8]);
        table[10] = mf.two_times(table[9]);

        table[20] = mf.square(table[10]);

        table[11] = mf.two_times(table[10]);
        table[12] = mf.two_times(table[11]);
        table[13] = mf.two_times(table[12]);
        table[14] = mf.two_times(table[13]);

        table[28] = mf.square(table[14]);
        table[15] = mf.two_times(table[14]);

        table[16] = mf.two_times(table[15]);
        table[21] = mf.two_times(table[20]);
        table[17] = mf.two_times(table[16]);
        table[22] = mf.two_times(table[21]);
        table[18] = mf.two_times(table[17]);
        table[23] = mf.two_times(table[22]);
        table[19] = mf.two_times(table[18]);
        table[24] = mf.two_times(table[23]);

        table[25] = mf.two_times(table[24]);
        table[29] = mf.two_times(table[28]);
        table[26] = mf.two_times(table[25]);
        table[30] = mf.two_times(table[29]);
        table[27] = mf.two_times(table[26]);
        table[31] = mf.two_times(table[30]);
  } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 8) {
        table[9]  = mf.two_times(table[8]);

        table[18] = mf.square(table[9]);
        table[10] = mf.two_times(table[9]);
        table[11] = mf.two_times(table[10]);
        table[12] = mf.two_times(table[11]);

        table[24] = mf.square(table[12]);
        table[13] = mf.two_times(table[12]);
        table[14] = mf.two_times(table[13]);

        table[29] = mf.multiply(table[18], table[11]);

        table[15] = mf.two_times(table[14]);
        table[19] = mf.two_times(table[18]);
        table[16] = mf.two_times(table[15]);
        table[20] = mf.two_times(table[19]);
        table[17] = mf.two_times(table[16]);
        table[21] = mf.two_times(table[20]);

        table[22] = mf.two_times(table[21]);
        table[25] = mf.two_times(table[24]);
        table[23] = mf.two_times(table[22]);
        table[26] = mf.two_times(table[25]);

        table[27] = mf.two_times(table[26]);
        table[30] = mf.two_times(table[29]);
        table[28] = mf.two_times(table[27]);
        table[31] = mf.two_times(table[30]);
  } else {
        table[9]  = mf.two_times(table[8]);
        table[10] = mf.two_times(table[9]);
        table[11] = mf.two_times(table[10]);
        table[12] = mf.two_times(table[11]);

        table[24] = mf.square(table[12]);
        table[13] = mf.two_times(table[12]);
        table[14] = mf.two_times(table[13]);
        table[15] = mf.two_times(table[14]);
        table[16] = mf.two_times(table[15]);

        table[17] = mf.two_times(table[16]);
        table[25] = mf.two_times(table[24]);

        table[18] = mf.two_times(table[17]);
        table[26] = mf.two_times(table[25]);

        table[19] = mf.two_times(table[18]);
        table[27] = mf.two_times(table[26]);

        table[20] = mf.two_times(table[19]);
        table[28] = mf.two_times(table[27]);

        table[21] = mf.two_times(table[20]);
        table[29] = mf.two_times(table[28]);

        table[22] = mf.two_times(table[21]);
        table[30] = mf.two_times(table[29]);

        table[23] = mf.two_times(table[22]);
        table[31] = mf.two_times(table[30]);
  }
}
    }
    else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 64) {

        // for MontyHalfRange and MontyFullRangeMasked, the two_times calls
        // would be slightly faster if table were type C (and using
        // getCanonicalValue as needed), but I don't expect these code sections
        // to ever be fast enough to be worth optimizing with table type C.
        table[0] = mf.getUnityValue();   // montgomery one
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        for (size_t i=0; i<63; ++i)
            table[i+1] = mf.two_times(table[i]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<23; ++i)
            table[i+1] = mf.two_times(table[i]);
        table[46] = mf.square(table[23]);
        table[24] = mf.two_times(table[23]);
        table[25] = mf.two_times(table[24]);
        table[26] = mf.two_times(table[25]);
        table[27] = mf.two_times(table[26]);
        table[28] = mf.two_times(table[27]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=28; i<45; ++i) {
            table[i+1] = mf.two_times(table[i]);
            table[i+19] = mf.two_times(table[i+18]);
        }
} else {
// This tested fastest for me when using __uint128_t
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<16; ++i)
            table[i+1] = mf.two_times(table[i]);

        table[32] = mf.square(table[16]);

        table[17] = mf.two_times(table[16]);
        table[18] = mf.two_times(table[17]);
        table[19] = mf.two_times(table[18]);
        table[20] = mf.two_times(table[19]);
        table[21] = mf.two_times(table[20]);

        table[22] = mf.two_times(table[21]);
        table[33] = mf.two_times(table[32]);
        table[53] = mf.multiply(table[21], table[32]);

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=22; i<31; ++i) {
            table[i+1] = mf.two_times(table[i]);
            table[i+12] = mf.two_times(table[i+11]);
        }
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=42; i<52; ++i) {
            table[i+1] = mf.two_times(table[i]);
            table[i+12] = mf.two_times(table[i+11]);
        }
}
    }
    else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 128) {

        // for MontyHalfRange and MontyFullRangeMasked, the two_times calls
        // would be slightly faster if table were type C (and using
        // getCanonicalValue as needed), but I don't expect these code sections
        // to ever be fast enough to be worth optimizing with table type C.
        table[0] = mf.getUnityValue();   // montgomery one
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        HPBC_CLOCKWORK_ASSERT2(TABLESIZE % 128 == 0);
        for (size_t i=0; i<TABLESIZE-1; ++i)
            table[i+1] = mf.two_times(table[i]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<15; ++i)
            table[i+1] = mf.two_times(table[i]);
        table[30] = mf.square(table[15]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=15; i<23; ++i)
            table[i+1] = mf.two_times(table[i]);
        // 0 -> 23, 30

        // precondition 30
        table[46] = mf.square(table[23]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=23; i<27; ++i) {
            table[i+1] = mf.two_times(table[i]);     // 24 -> 27
            table[i+8] = mf.two_times(table[i+7]);   // 31 -> 34
        }
        // 0 -> 27, 30 -> 34, 46

        // precondition 46, 34
        table[68] = mf.square(table[34]);
        table[28] = mf.two_times(table[27]);
        table[35] = mf.two_times(table[34]);
        table[29] = mf.two_times(table[28]);
        table[36] = mf.two_times(table[35]);
        // 0 -> 36, 46, 68

        for (size_t i=36; i<45; ++i) {
            table[i+1] = mf.two_times(table[i]);     // 37 -> 45
            table[i+11] = mf.two_times(table[i+10]); // 47 -> 55
        }
        // 0 -> 55, 68

        // precondition 55, 68, 48
        table[56] = mf.two_times(table[55]);
        table[69] = mf.two_times(table[68]);
        // 0 -> 56, 68 -> 69
        table[96] = mf.square(table[48]);
        table[112] = mf.square(table[56]);
        for (size_t i=56; i<67; ++i) {
            table[i+1] = mf.two_times(table[i]);      // 57 -> 67
            table[i+14] = mf.two_times(table[i+13]);  // 70 -> 80
        }
        // 0 -> 80, 96, 112

        // needs 80 96 112
        for (size_t i=80; i<95; ++i) {
            table[i+1] = mf.two_times(table[i]);      // 81 -> 95
            table[i+17] = mf.two_times(table[i+16]);  // 97 -> 111
            table[i+33] = mf.two_times(table[i+32]);  // 113 -> 127
        }
} else {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<12; ++i)
            table[i+1] = mf.two_times(table[i]);
        // we now have 0 -> 12 (inclusive)

        table[24] = mf.square(table[12]);

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=12; i<18; ++i)
            table[i+1] = mf.two_times(table[i]);
        // 0 -> 18, and 24

        table[36] = mf.square(table[18]);
        table[42] = mf.multiply(table[18], table[24]);

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=18; i<23; ++i) {
            table[i+1] = mf.two_times(table[i]);
            table[i+7] = mf.two_times(table[i+6]);
        }
        // 0 -> 29, and 36 and 42

        table[72] = mf.square(table[36]);
        table[84] = mf.square(table[42]);
        
        for (size_t i=29; i<34; ++i) {
            table[i+1] = mf.two_times(table[i]);
            table[i+8] = mf.two_times(table[i+7]);
            table[i+14] = mf.two_times(table[i+13]);
        }
        // 0 -> 34, 36 -> 47, 72, 84
        table[35] = mf.two_times(table[34]);
        table[48] = mf.two_times(table[47]);
        // 0 -> 48, 72, 84

        table[106] = mf.multiply(table[34], table[72]);
        table[117] = mf.multiply(table[45], table[72]);

        for (size_t i=48; i<59; ++i) {
            table[i+1] = mf.two_times(table[i]);
            table[i+25] = mf.two_times(table[i+24]);
            table[i+37] = mf.two_times(table[i+36]);
        }
        // 0->59, 72->95, 106, 117

        for (size_t i=59; i<69; ++i) {
            table[i+1] = mf.two_times(table[i]);     // 60 -> 69
            table[i+37] = mf.two_times(table[i+36]); // 96 -> 105
            table[i+48] = mf.two_times(table[i+47]); // 107 -> 116
            table[i+59] = mf.two_times(table[i+58]); // 118 -> 127
        }
        // 0 -> 69, 72 -> 105, 106, 107 -> 116, 117, 118 -> 127
        table[70] = mf.two_times(table[69]);
        table[71] = mf.two_times(table[70]);
}
    }
    else {
        // for MontyHalfRange and MontyFullRangeMasked, the two_times calls
        // would be slightly faster if table were type C (and using
        // getCanonicalValue if needed), but I don't expect this code clause to
        // ever be fast enough to be worth optimizing with table type C.

        HPBC_CLOCKWORK_ASSERT2(TABLESIZE % 256 == 0);
        table[0] = mf.getUnityValue();   // montgomery one
        for (size_t i=0; i<TABLESIZE-1; ++i)
            table[i+1] = mf.two_times(table[i]);
    }



    constexpr size_t MASK = TABLESIZE - 1;
    using RU = typename MFE::RU;
    constexpr int digitsRU = hc::ut_numeric_limits<RU>::digits;
    constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));

    // recall that we set  constexpr int P = static_cast<int>(TABLE_BITS);

    int shift;
    V result;
    if HURCHALLA_CPP17_CONSTEXPR (P2 > P) {
        constexpr size_t MASKBIG = (1u << P2) - 1u;
        if (n <= MASKBIG) {
            size_t loindex = static_cast<size_t>(n);
            HPBC_CLOCKWORK_ASSERT2(loindex < ut_numeric_limits<RU>::digits);
            result = MFE::twoPowLimited(mf, loindex);
            return result;
        }

        // count_leading_zeros returns the number of leading 0-bits in n, starting
        // at the most significant bit position. If n is 0, the result is undefined
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;

        // because we returned above if (n <= MASKBIG), we can assert the following:
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);
        shift = numbits - P2;

        U tmp = n >> shift;
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        // normally we'd use (tmp & MASKBIG), but it's redundant with tmp <= MASKBIG
        size_t index = static_cast<size_t>(tmp);
        HPBC_CLOCKWORK_ASSERT2(index < ut_numeric_limits<RU>::digits);
        result = MFE::twoPowLimited(mf, index);
    }
    else {
        if (n <= MASK) {
            result = table[static_cast<size_t>(n)];
            return result;
        }

        // count_leading_zeros returns the number of leading 0-bits in n, starting
        // at the most significant bit position. If n is 0, the result is undefined
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        // because we returned above if (n <= MASK), we can assert the following:
        HPBC_CLOCKWORK_ASSERT2(numbits > P);

        shift = numbits - P;
        U tmp = n >> shift;
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
        // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
        size_t index = static_cast<size_t>(tmp);
        result = table[index];
    }
    HPBC_CLOCKWORK_ASSERT2(shift > 0);


    while (shift >= P) {
        if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
            while (shift > P && (static_cast<size_t>(n>>(shift-1)) & 1) == 0) {
                result = mf.square(result);
                --shift;
            }
        }

        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            SV sv = MFE::getSquaringValue(mf, result);
            HPBC_CLOCKWORK_ASSERT2(P > 0);
            HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P - 1; ++i)
                sv = MFE::squareSV(mf, sv);
            result = MFE::squareToMontgomeryValue(mf, sv);
        } else {
            HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P; ++i)
                result = mf.square(result);
        }

        shift -= P;
        size_t index = static_cast<size_t>(n >> shift) & MASK;
        result = mf.multiply(result, table[index]);
    }


    if (shift == 0)
        return result;
    HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P);

    for (int i=0; i<shift; ++i)
        result = mf.square(result);
    size_t tmpmask = (1u << shift) - 1;
    size_t index = static_cast<size_t>(n) & tmpmask;
    result = mf.multiply(result, table[index]);
    return result;
  }





#ifdef __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpass-failed"
#endif
  
  // Array version of montgomery two pow
  template <class MF, typename U,
            size_t ARRAY_SIZE, size_t TABLE_BITS = 5, size_t CODE_SECTION = 0,
            bool USE_SQUARING_VALUE_OPTIMIZATION = false>
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    // FYI you almost certainly want either TABLE_BITS == 0, or
    // 2 <= TABLE_BITS <= 6.  Anything above 9 is probably a very bad idea even
    // if it works since it would mean calculating 1024+ table entries!
    static_assert(0 <= TABLE_BITS && TABLE_BITS < 10, "");

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
        if (n_max < n[j])
            n_max = n[j];
    }


    if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        // an array version of scalar two_pow TABLE_BITS 0, CODE_SECTION 2

        std::array<V, ARRAY_SIZE> result;
        if (n_max <= MASK) {
            for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j]= MFE_LU::twoPowLimited(mf[j], static_cast<size_t>(n[j]));
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);

        int shift = numbits - P2;
        std::array<size_t, ARRAY_SIZE> tmp;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            tmp[j] = static_cast<size_t>(n[j] >> shift);
            HPBC_CLOCKWORK_ASSERT2(tmp[j] <= MASK);
            // normally we use (tmp & MASK), but it's redundant with tmp <= MASK
            result[j] = MFE_LU::twoPowLimited(mf[j], tmp[j]);
        }

        while (shift >= P2) {
            shift -= P2;
            std::array<size_t, ARRAY_SIZE> index;
            std::array<V, ARRAY_SIZE> tableVal;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                tmp[j] = static_cast<size_t>(n[j] >> shift);
                index[j] = tmp[j] & MASK;
                tableVal[j] = MFE_LU::twoPowLimited(mf[j], index[j]);
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

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   tableVal[j]);
            }
        }
        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);
        size_t tmpmask = (1u << shift) - 1u;

        std::array<size_t, ARRAY_SIZE> index;
        std::array<V, ARRAY_SIZE> tableVal;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            index[j] = static_cast<size_t>(n[j]) & tmpmask;
            tableVal[j] = MFE_LU::twoPowLimited(mf[j], index[j]);
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

    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
        // an array version of scalar two_pow TABLE_BITS 0, CODE_SECTION 17
        // (CODE_SECTION 17 is an updated version of CODE_SECTION 3)

        constexpr int NUMBITS_MASKBIG = P2 + 1;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

        std::array<C, ARRAY_SIZE> cR1;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            cR1[j] = MFE_LU::getMontvalueR(mf[j]);

        std::array<V, ARRAY_SIZE> result;
        if (n_max <= MASK) {
            for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = MFE_LU::twoPowLimited_times_x(mf[j], static_cast<size_t>(n[j]), cR1[j]);
            return result;
        }

        std::array<C, ARRAY_SIZE> cR2;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            cR2[j] = mf[j].getCanonicalValue(mf[j].template square<LowuopsTag>(cR1[j]));

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits >= NUMBITS_MASKBIG);

        int shift = numbits - NUMBITS_MASKBIG;
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t tmp = static_cast<size_t>(n[j] >> shift);
            HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
            size_t loindex = tmp & MASK;
            C cHigh = cR1[j];
            cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR2[j]);
            result[j] = MFE_LU::twoPowLimited_times_x(mf[j], loindex, cHigh);
        }

        while (shift >= NUMBITS_MASKBIG) {
            shift -= NUMBITS_MASKBIG;

            std::array<V, ARRAY_SIZE> val1;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t tmp = static_cast<size_t>(n[j] >> shift);
                size_t loindex = tmp & MASK;
                C cHigh = cR1[j];
                cHigh.cmov(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u), cR2[j]);
                val1[j] = MFE_LU::twoPowLimited_times_x(mf[j], loindex, cHigh);
            }

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                std::array<SV, ARRAY_SIZE> sv;
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::getSquaringValue(mf[j], result[j]);

                static_assert(NUMBITS_MASKBIG > 0, "");
                for (int i=0; i<NUMBITS_MASKBIG - 1; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf[j], sv[j]);
                }
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = MFE_LU::squareToMontgomeryValue(mf[j], sv[j]);
            } else {
                for (int i=0; i<NUMBITS_MASKBIG; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                }
            }

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                    val1[j]);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);
        size_t tmpmask = (1u << shift) - 1u;

        std::array<V, ARRAY_SIZE> val1;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t tmp = static_cast<size_t>(n[j]) & tmpmask;
            HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
            HPBC_CLOCKWORK_ASSERT2(((tmp >> (NUMBITS_MASKBIG - 1)) & 1u) == 0u);
            size_t loindex = tmp;
            C cHigh = cR1[j];
            val1[j] = MFE_LU::twoPowLimited_times_x(mf[j], loindex, cHigh);
        }

        for (int i=0; i<shift; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   val1[j]);
        }
        return result;

    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2 || CODE_SECTION == 3) {
        // an array version of scalar two_pow TABLE_BITS 1, CODE_SECTION 3 and 4
        // it's a simplified form of it too.
        int shift;
        std::array<V, ARRAY_SIZE> result;
        {
            if (n_max <= MASK) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                    size_t loindex = static_cast<size_t>(n[j]);
                    result[j] = MFE_LU::twoPowLimited(mf[j], loindex);
                }
                return result;
            }

            HPBC_CLOCKWORK_ASSERT2(n_max > 0);
            int leading_zeros = count_leading_zeros(n_max);
            int numbits = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > P2);

            shift = numbits - P2;
            HPBC_CLOCKWORK_ASSERT2(shift > 0);

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t tmp = static_cast<size_t>(n[j] >> shift);
                HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
                result[j] = MFE_LU::twoPowLimited(mf[j], tmp);
            }
        }

        while (shift >= 1) {
            --shift;
            if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2) {
                // branch based code
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                    result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                    if (static_cast<size_t>(n[j] >> shift) & 1u)
                        result[j] = mf[j].two_times(result[j]);
                }
            } else {
                // cmov based code
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                    result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                    V vtmp = mf[j].two_times(result[j]);
                       // result[j] = ((n[j] >> shift) & 1u) ? vtmp : result[j];
                    result[j].cmov(static_cast<size_t>(n[j] >> shift) & 1u, vtmp);
                }
            }
        }
        return result;

    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION >= 4 && CODE_SECTION <= 7) {
        // this is an array version of scalar two_pow TABLE_BITS 0,
        // CODE_SECTION 10-13

        // make the high table size either 2, 4, 8, or 16, depending on the CODE_SECTION
        constexpr int NUMBITS_TABLE_HIGH_SIZE = static_cast<int>(CODE_SECTION) - 3;

        static_assert(1 <= NUMBITS_TABLE_HIGH_SIZE && NUMBITS_TABLE_HIGH_SIZE <= 4, "");

        constexpr size_t TABLE_HIGH_SIZE = 1u << NUMBITS_TABLE_HIGH_SIZE;
        constexpr int NUMBITS_MASKBIG = P2 + NUMBITS_TABLE_HIGH_SIZE;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

        std::array<std::array<C, ARRAY_SIZE>, TABLE_HIGH_SIZE> table_high;


        static_assert(2 <= TABLE_HIGH_SIZE && TABLE_HIGH_SIZE <= 16, "");

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            C cR1 = MFE_LU::getMontvalueR(mf[j]);
            V vR2 = mf[j].template square<hc::LowuopsTag>(cR1);
            table_high[0][j] = cR1;                                // R^1
            table_high[1][j] = mf[j].getCanonicalValue(vR2);       // R^2
        }

        if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE >= 4) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR3 = mf[j].template multiply<hc::LowuopsTag>(table_high[1][j], table_high[0][j]);
                table_high[2][j] = mf[j].getCanonicalValue(vR3);   // R^3
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR4 = mf[j].template square<hc::LowuopsTag>(table_high[1][j]);
                table_high[3][j] = mf[j].getCanonicalValue(vR4);   // R^4
            }
        }

        if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE >= 8) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR5 = mf[j].template multiply<hc::LowuopsTag>(table_high[2][j], table_high[1][j]);
                table_high[4][j] = mf[j].getCanonicalValue(vR5);   // R^5
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR6 = mf[j].template square<hc::LowuopsTag>(table_high[2][j]);
                table_high[5][j] = mf[j].getCanonicalValue(vR6);   // R^6
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR7 = mf[j].template multiply<hc::LowuopsTag>(table_high[3][j], table_high[2][j]);
                table_high[6][j] = mf[j].getCanonicalValue(vR7);   // R^7
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR8 = mf[j].template square<hc::LowuopsTag>(table_high[3][j]);
                table_high[7][j] = mf[j].getCanonicalValue(vR8);   // R^8
            }
        }

        if HURCHALLA_CPP17_CONSTEXPR (TABLE_HIGH_SIZE >= 16) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR9 = mf[j].template multiply<hc::LowuopsTag>(table_high[4][j], table_high[3][j]);
                table_high[8][j] = mf[j].getCanonicalValue(vR9);     // R^9
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR10 = mf[j].template square<hc::LowuopsTag>(table_high[4][j]);
                table_high[9][j] = mf[j].getCanonicalValue(vR10);    // R^10
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR11 = mf[j].template multiply<hc::LowuopsTag>(table_high[5][j], table_high[4][j]);
                table_high[10][j] = mf[j].getCanonicalValue(vR11);   // R^11
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR12 = mf[j].template square<hc::LowuopsTag>(table_high[5][j]);
                table_high[11][j] = mf[j].getCanonicalValue(vR12);   // R^12
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR13 = mf[j].template multiply<hc::LowuopsTag>(table_high[6][j], table_high[5][j]);
                table_high[12][j] = mf[j].getCanonicalValue(vR13);   // R^13
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR14 = mf[j].template square<hc::LowuopsTag>(table_high[6][j]);
                table_high[13][j] = mf[j].getCanonicalValue(vR14);   // R^14
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR15 = mf[j].template multiply<hc::LowuopsTag>(table_high[7][j], table_high[6][j]);
                table_high[14][j] = mf[j].getCanonicalValue(vR15);   // R^15
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                V vR16 = mf[j].template square<hc::LowuopsTag>(table_high[7][j]);
                table_high[15][j] = mf[j].getCanonicalValue(vR16);   // R^16
            }
        }

        int shift = 0;
        if (n_max > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n_max > 0);
            int leading_zeros = count_leading_zeros(n_max);
            int numbits = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
        }
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);

        std::array<V, ARRAY_SIZE> result;

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t tmp = static_cast<size_t>(n[j] >> shift);
            HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
            size_t loindex = tmp & MASK;
            size_t hiindex = tmp >> P2;
            HPBC_CLOCKWORK_ASSERT2(hiindex < TABLE_HIGH_SIZE);
            result[j] = MFE_LU::twoPowLimited_times_x(mf[j], loindex, table_high[hiindex][j]);
        }

        while (shift >= NUMBITS_MASKBIG) {
            shift -= NUMBITS_MASKBIG;

            std::array<V, ARRAY_SIZE> val1;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t tmp = static_cast<size_t>(n[j] >> shift);
                size_t loindex = tmp & MASK;
                size_t hiindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
                val1[j] = MFE_LU::twoPowLimited_times_x(mf[j], loindex, table_high[hiindex][j]);
            }

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                std::array<SV, ARRAY_SIZE> sv;
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::getSquaringValue(mf[j], result[j]);

                static_assert(NUMBITS_MASKBIG > 0, "");
                for (int i=0; i<NUMBITS_MASKBIG - 1; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf[j], sv[j]);
                }
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = MFE_LU::squareToMontgomeryValue(mf[j], sv[j]);
            } else {
                for (int i=0; i<NUMBITS_MASKBIG; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                }
            }

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                    val1[j]);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);
        size_t tmpmask = (1u << shift) - 1u;

        std::array<V, ARRAY_SIZE> val1;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t tmp = static_cast<size_t>(n[j]) & tmpmask;
            size_t loindex = tmp & MASK;
            size_t hiindex = (tmp >> P2) & (TABLE_HIGH_SIZE - 1);
            val1[j] = MFE_LU::twoPowLimited_times_x(mf[j], loindex, table_high[hiindex][j]);
        }

        for (int i=0; i<shift; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   val1[j]);
        }
        return result;

    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 8) {
        // standard k-ary pow with table init

        constexpr int P = static_cast<int>(TABLE_BITS);
        static_assert(P >= 0, "");
        constexpr size_t TABLESIZE = 1u << P;   
        static_assert(TABLESIZE >= 1, "");
        constexpr size_t TABLE_MASK = TABLESIZE - 1u;

        int shift;
        std::array<V, ARRAY_SIZE> result;
        {
            if (n_max <= MASK) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                    size_t loindex = static_cast<size_t>(n[j]);
                    result[j] = MFE_LU::twoPowLimited(mf[j], loindex);
                }
                return result;
            }

            HPBC_CLOCKWORK_ASSERT2(n_max > 0);
            int leading_zeros = count_leading_zeros(n_max);
            int numbits = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > P2);

            shift = numbits - P2;
            HPBC_CLOCKWORK_ASSERT2(shift > 0);

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t tmp = static_cast<size_t>(n[j] >> shift);
                HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
                result[j] = MFE_LU::twoPowLimited(mf[j], tmp);
            }
        }

        // Initialize the precalculation table...
        // We'll likely get very good instruction level parallelism via the array,
        // and have no need for any of the init tricks seen in the scalar two_pow 
        // functions above.
        C table[TABLESIZE][ARRAY_SIZE];
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[0][j] = mf[j].getUnityValue();   // montgomery one
        for (size_t i=1; i<TABLESIZE; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i][j] = mf[j].two_times(table[i-1][j]);
        }

        while (shift >= P) {
            shift -= P;

            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                std::array<SV, ARRAY_SIZE> sv;
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::getSquaringValue(mf[j], result[j]);
                static_assert(P > 0, "");
                for (int i=0; i<P - 1; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf[j], sv[j]);
                }
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = MFE_LU::squareToMontgomeryValue(mf[j], sv[j]);
            } else {
                for (int i=0; i<P; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                }
            }

            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t tmp = static_cast<size_t>(n[j] >> shift);
                size_t index = tmp & TABLE_MASK;
                result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   table[index][j]);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P);

        for (int i=0; i<shift; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }
        size_t tmpmask = (1u << shift) - 1u;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n[j]) & tmpmask;
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                            table[index][j]);
        }
        return result;

    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 27) {
        // this corresponds to scalar two_pow's CODE_SECTION 27

        std::array<V, ARRAY_SIZE> result;
        if (n_max <= MASK) {
            for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j]= MFE_LU::twoPowLimited(mf[j], static_cast<size_t>(n[j]));
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);

        int shift = numbits - P2;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n[j] >> shift);
            HPBC_CLOCKWORK_ASSERT2(index <= MASK);
            // normally we use (index & MASK), but it's redundant with index <= MASK
            result[j] = MFE_LU::twoPowLimited(mf[j], index);
        }

        while (shift >= P2) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                result[j] = mf[j].two_times(result[j]);
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
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                size_t index = static_cast<size_t>(n[j] >> shift) & MASK;
                C tmp = mf[j].getCanonicalValue(result[j]);
                result[j] = MFE_LU::twoPowLimited_times_x(mf[j], index, tmp);
            }
        }
        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P2);
        size_t tmpmask = (1u << shift) - 1u;

        std::array<V, ARRAY_SIZE> tableVal;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n[j]) & tmpmask;
            tableVal[j] = MFE_LU::twoPowLimited(mf[j], index);
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
    } else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 28) {
        // this corresponds to scalar two_pow's CODE_SECTION 28

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
            size_t index = static_cast<size_t>(n[j] >> shift);
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
                size_t index = static_cast<size_t>(n[j] >> shift) & MASK;
                C tmp = mf[j].getCanonicalValue(result[j]);
                result[j] = MFE_LU::twoPowLimited_times_x_v2(mf[j], index + 1, tmp);
            }
        }

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            C tmp = mf[j].getCanonicalValue(result[j]);
            result[j] = mf[j].template divideBySmallPowerOf2<hc::LowuopsTag>(tmp, 1);
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
    } else {
        // this corresponds to scalar two_pow's CODE_SECTION 29
        static_assert(CODE_SECTION == 29, "");

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
                size_t index = static_cast<size_t>(n[j] >> shift) & MASK;
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
            size_t index = static_cast<size_t>(n[j] >> shift) & MASK;
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
    }

  }


#ifdef __clang__
#  pragma GCC diagnostic pop
#endif
};


}} // end namespace

#endif
