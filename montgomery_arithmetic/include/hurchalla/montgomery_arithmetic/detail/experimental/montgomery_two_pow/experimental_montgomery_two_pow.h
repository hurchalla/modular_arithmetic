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
            size_t TABLE_BITS = 5, size_t CODE_SECTION = 0>
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
        if (n <= MASK) {
            size_t index = static_cast<size_t>(n);
            RU num = static_cast<RU>(static_cast<RU>(1) << index);
            return MFE::convertInExtended(mf, num);
        }

        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);

        int shift = numbits - P2;
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
            HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2; ++i)
                result = mf.square(result);
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
        if (n <= MASK) {
            size_t loindex = static_cast<size_t>(n);
            RU num = static_cast<RU>(static_cast<RU>(1) << loindex);
            return MFE::convertInExtended(mf, num);
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
        V result = MFE::convertInExtended_aTimesR(mf, num, magicValue);

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

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                    result = mf.square(result);

                result = mf.multiply(result, val1);
            }
            else {
                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                num = static_cast<RU>(static_cast<RU>(1) << loindex);
                V val1 = MFE::convertInExtended_aTimesR(mf, num, magicValue);
                V val2 = MFE::convertInExtended(mf, num);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                    result = mf.square(result);

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
        if (n <= MASK)
            return MFE::twoPowLimited(mf, static_cast<size_t>(n));
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);

        int shift = numbits - P2;
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
            HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P2; ++i)
                result = mf.square(result);
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
} else {
        // This is basically a copy of code section 1, except we replace calls
        // to convertInExtended() with twoPowLimited(), and we replace calls to
        // convertInExtended_aTimesR() with RTimesTwoPowLimited().
        if (n <= MASK) {
            size_t loindex = static_cast<size_t>(n);
            return MFE::twoPowLimited(mf, loindex);
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
        V result = MFE::RTimesTwoPowLimited(mf, loindex, magicValue);

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

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                    result = mf.square(result);

                result = mf.multiply(result, val1);
            }
            else {
                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                V val1 = MFE::RTimesTwoPowLimited(mf, loindex, magicValue);
                V val2 = MFE::twoPowLimited(mf, loindex);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<(P2 + 1); ++i)
                    result = mf.square(result);

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
}
    }
    else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 2) {
        table[0] = mf.getUnityValue();   // montgomery one

// The different code sections should be functionally equivalent.  You can test
// to see which is fastest.
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        table[1] = mf.two_times(table[0]);
} else {
// For this particular case of TABLE_BITS == 1 (TABLESIZE == 2), we can
// use a version of 2^k-ary that is heavily optimized for the 1 bit table:
        if (n <= 1) {
            if (n == 0)
                return mf.getUnityValue();
            else
                return mf.two_times(mf.getUnityValue());
        }
        HPBC_CLOCKWORK_ASSERT2(n > 1);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > 1);
        int shift = numbits - 1;
        HPBC_CLOCKWORK_ASSERT2((n >> shift) == 1);   // because shift == numbits - 1
        C mont_two = mf.two_times(mf.getUnityValue());
        C cresult = mont_two;
        HPBC_CLOCKWORK_ASSERT2(shift >= 1);
        // since cresult == 2, this two_times() is equivalent to squaring
        cresult = mf.two_times(cresult);
        --shift;
        V result;
  if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
    // branch based code
        if (static_cast<size_t>(n >> shift) & 1)
            cresult = mf.two_times(cresult);
        result = cresult;
        while (shift >= 1) {
            result = mf.square(result);
            --shift;
            if (static_cast<size_t>(n >> shift) & 1)
                result = mf.two_times(result);
        }
  } else {
    // cmov based code
        C ctmp = mf.two_times(cresult);
        cresult.cmov(static_cast<size_t>(n >> shift) & 1, ctmp);
        result = cresult;
        while (shift >= 1) {
            result = mf.square(result);
            --shift;
            V vtmp = mf.two_times(result);
            result.cmov(static_cast<size_t>(n >> shift) & 1, vtmp);
        }
  }
        return result;
}
    } else if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 4) {
        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
} else {
// try the same optimization as TABLE_BITS == 1 (TABLESIZE == 2), but
// using the 4 entry table that we have here.

        // this section is a copy/paste of the main code at bottom of this
        // function.  Note P == 2 here.
        //
        constexpr size_t MASK = TABLESIZE - 1;
        if (n <= MASK)
            return table[static_cast<size_t>(n)];
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P);
        int shift = numbits - P;
        U tmp = n >> shift;
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
        size_t index = static_cast<size_t>(tmp);
        V result = table[index];

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
        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = mf.two_times(table[0]);
        table[2] = mf.two_times(table[1]);
        table[3] = mf.two_times(table[2]);
        table[4] = mf.two_times(table[3]);
        table[5] = mf.two_times(table[4]);
        table[6] = mf.two_times(table[5]);
        table[7] = mf.two_times(table[6]);
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
} else {
// try the same optimization as TABLE_BITS == 1 (TABLESIZE == 2), but
// using the 8 entry table that we have here.

        // this section is a copy/paste of the main code at bottom of this
        // function.  Note P == 3 here.
        //
        constexpr size_t MASK = TABLESIZE - 1;
        if (n <= MASK)
            return table[static_cast<size_t>(n)];
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P);
        int shift = numbits - P;
        U tmp = n >> shift;
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
        size_t index = static_cast<size_t>(tmp);
        V result = table[index];

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
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        table[0] = mf.getUnityValue();   // montgomery one
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<15; ++i)
            table[i+1] = mf.two_times(table[i]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
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
if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        table[0] = mf.getUnityValue();   // montgomery one
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<31; ++i)
            table[i+1] = mf.two_times(table[i]);
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
//        note: if using multiply, it would have to be moved lower
//        table[19] = mf.multiply(table[9]), table[10]);
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
        table[16] = MFE::twoPowLimited(mf, 16);
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
    // convertIn(), but it's dubious to do so - it ties us to using a
    // MontgomeryForm that computes RSquaredModN in its constructor.
    // Thus we don't call convertIn, even though (if we already had computed
    // RSquaredModN) it might in some cases be faster.
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

        table[29] = mf.multiply(table[14], table[15]);

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
        table[53] = mf.multiply(table[32], table[21]);

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
        table[0] = mf.getUnityValue();   // montgomery one
        //return mf.getUnityValue();
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
        table[42] = mf.multiply(table[24], table[18]);

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

        table[106] = mf.multiply(table[72], table[34]);
        table[117] = mf.multiply(table[72], table[45]);

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
        HPBC_CLOCKWORK_ASSERT2(TABLESIZE % 256 == 0);
        table[0] = mf.getUnityValue();   // montgomery one
        for (size_t i=0; i<TABLESIZE-1; ++i)
            table[i+1] = mf.two_times(table[i]);
    }


    constexpr size_t MASK = TABLESIZE - 1;
    if (n <= MASK)
        return table[static_cast<size_t>(n)];


    // count_leading_zeros returns the number of leading 0-bits in n, starting
    // at the most significant bit position. If n is 0, the result is undefined
    HPBC_CLOCKWORK_ASSERT2(n > 0);
    int leading_zeros = count_leading_zeros(n);
    int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
    // because we returned above if (n <= MASK), we can assert the following:
    HPBC_CLOCKWORK_ASSERT2(numbits > P);

    int shift = numbits - P;
    U tmp = n >> shift;
    HPBC_CLOCKWORK_ASSERT2(tmp <= MASK);
    // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
    size_t index = static_cast<size_t>(tmp);
    V result = table[index];


    while (shift >= P) {
        if (USE_SLIDING_WINDOW_OPTIMIZATION) {
            while (shift > P && (static_cast<size_t>(n>>(shift-1)) & 1) == 0) {
                result = mf.square(result);
                --shift;
            }
        }

//        HURCHALLA_REQUEST_UNROLL_LOOP
        for (int i=0; i<P; ++i)
            result = mf.square(result);

        shift -= P;
        index = static_cast<size_t>(n >> shift) & MASK;
        result = mf.multiply(result, table[index]);
    }


    if (shift == 0)
        return result;
    HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < P);

    for (int i=0; i<shift; ++i)
        result = mf.square(result);
    size_t tmpmask = (1u << shift) - 1;
    index = static_cast<size_t>(n) & tmpmask;
    result = mf.multiply(result, table[index]);
    return result;
  }





#ifdef __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpass-failed"
#endif
  
  // Array version of montgomery two pow
  template <class MF, typename U,
            size_t ARRAY_SIZE, size_t TABLE_BITS = 5, size_t CODE_SECTION = 0>
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
    using MFE = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;
    using std::size_t;

    constexpr int P = static_cast<int>(TABLE_BITS);

    // Initialize the precalculation table...
    // We'll likely get very good instruction level parallelism via the array,
    // and have no need for any of the init tricks seen in the function above.
    static_assert(P >= 0, "");
    constexpr size_t TABLESIZE = 1u << P;   
    static_assert(TABLESIZE >= 1, "");

    U n_max = n[0];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=1; j<ARRAY_SIZE; ++j) {
        if (n_max < n[j])
            n_max = n[j];
    }


    if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE == 1) {

if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {

        using RU = typename MFE::RU;
        constexpr int digitsRU = hc::ut_numeric_limits<RU>::digits;
        constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));
        constexpr size_t MASK = (1u << P2) - 1u;

        if (n_max <= MASK) {
            std::array<V, ARRAY_SIZE> result;
            for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j]= MFE::twoPowLimited(mf[j], static_cast<size_t>(n[j]));
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > P2);

        int shift = numbits - P2;
        std::array<V, ARRAY_SIZE> result;
        std::array<size_t, ARRAY_SIZE> tmp;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            tmp[j] = static_cast<size_t>(n[j] >> shift);
            HPBC_CLOCKWORK_ASSERT2(tmp[j] <= MASK);
            // normally we use (tmp & MASK), but it's redundant with tmp <= MASK
            result[j] = MFE::twoPowLimited(mf[j], tmp[j]);
        }

        while (shift >= P2) {
            shift -= P2;
            std::array<size_t, ARRAY_SIZE> index;
            std::array<V, ARRAY_SIZE> tableVal;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                tmp[j] = static_cast<size_t>(n[j] >> shift);
                index[j] = tmp[j] & MASK;
                tableVal[j] = MFE::twoPowLimited(mf[j], index[j]);
            }

            for (int i=0; i<P2; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP for(size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j]= mf[j].template square<hc::LowuopsTag>(result[j]);
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
            tableVal[j] = MFE::twoPowLimited(mf[j], index[j]);
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

        using RU = typename MFE::RU;
        constexpr int digitsRU = hc::ut_numeric_limits<RU>::digits;
        constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));
        constexpr size_t MASK = (1u << P2) - 1u;

        if (n_max <= MASK) {
            std::array<V, ARRAY_SIZE> result;
            for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j]= MFE::twoPowLimited(mf[j], static_cast<size_t>(n[j]));
            return result;
        }
        std::array<RU, ARRAY_SIZE> magicValue;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            magicValue[j] = MFE::getMagicValue(mf[j]);
        }

        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits >= (P2 + 1));

        int shift = numbits - (P2 + 1);
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        std::array<V, ARRAY_SIZE> result;
        std::array<size_t, ARRAY_SIZE> tmp;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            tmp[j] = static_cast<size_t>(n[j] >> shift);
            HPBC_CLOCKWORK_ASSERT2(tmp[j] <= 2u*MASK + 1u);
            size_t loindex = tmp[j] & MASK;
            V val1 = MFE::RTimesTwoPowLimited(mf[j], loindex, magicValue[j]);
            V val2 = MFE::twoPowLimited(mf[j], loindex);
            size_t hibit = tmp[j] >> P2;
            HPBC_CLOCKWORK_ASSERT2(hibit <= 1);
            // val1 = (hibit == 0) ? val2 : val1;
            val1.cmov(hibit == 0, val2);
            result[j] = val1;
        }

        while (shift >= (P2 + 1)) {
            shift -= (P2 + 1);

            std::array<size_t, ARRAY_SIZE> loindex;
            std::array<V, ARRAY_SIZE> tableVal;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                tmp[j] = static_cast<size_t>(n[j] >> shift);
                loindex[j] = tmp[j] & MASK;
                V val1 = MFE::RTimesTwoPowLimited(mf[j], loindex[j],
                                                                 magicValue[j]);
                V val2 = MFE::twoPowLimited(mf[j], loindex[j]);
                size_t hibit = (tmp[j] >> P2) & 1u;
                //val1 = (hibit == 0) ? val2 : val1;
                val1.cmov(hibit == 0, val2);
                tableVal[j] = val1;
            }

            for (int i=0; i<(P2 + 1); ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP for(size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j]= mf[j].template square<hc::LowuopsTag>(result[j]);
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   tableVal[j]);
            }
        }
        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < (P2 + 1));
        size_t tmpmask = (1u << shift) - 1u;

        std::array<size_t, ARRAY_SIZE> index;
        std::array<V, ARRAY_SIZE> tableVal;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            index[j] = static_cast<size_t>(n[j]) & tmpmask;
            tableVal[j] = MFE::twoPowLimited(mf[j], index[j]);
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
else {
        C table[8][ARRAY_SIZE];
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[0][j] = mf[j].getUnityValue();   // montgomery one
        for (size_t i=1; i<8; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i][j] = mf[j].two_times(table[i-1][j]);
        }
        constexpr size_t MASK = 7;
        if (n_max <= MASK) {
            std::array<V, ARRAY_SIZE> result;
            for (size_t j=0; j<ARRAY_SIZE; ++j) {
                result[j] = table[static_cast<size_t>(n[j])][j];
            }
            return result;
        }
        HPBC_CLOCKWORK_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > 3);
        int shift = numbits - 3;

        std::array<V, ARRAY_SIZE> result;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n[j] >> shift);
            HPBC_CLOCKWORK_ASSERT2(index <= MASK);
            result[j] = table[index][j];
        }

        while (shift >= 1) {
            --shift;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
                V vtmp = mf[j].two_times(result[j]);
                   // result[j] = ((n[j] >> shift) & 1u) ? vtmp : result[j];
                result[j].cmov(static_cast<size_t>(n[j] >> shift) & 1u, vtmp);
            }
        }
        return result;
}
    }



    C table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
        table[0][j] = mf[j].getUnityValue();   // montgomery one
    for (size_t i=1; i<TABLESIZE; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            table[i][j] = mf[j].two_times(table[i-1][j]);
    }

    constexpr size_t MASK = TABLESIZE - 1u;
    if (n_max <= MASK) {
        std::array<V, ARRAY_SIZE> result;
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = table[static_cast<size_t>(n[j])][j];
        return result;
    }

    // count_leading_zeros returns the number of leading 0-bits in n_max,
    // starting at the most significant bit position. If n_max is 0, the
    // result is undefined
    HPBC_CLOCKWORK_ASSERT2(n_max > 0);
    int leading_zeros = count_leading_zeros(n_max);
    int numbits = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
    // because we returned above if (n_max <= MASK), we know the following:
    HPBC_CLOCKWORK_ASSERT2(numbits > P);

    int shift = numbits - P;
    std::array<V, ARRAY_SIZE> result;
    std::array<U, ARRAY_SIZE> tmp;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        tmp[j] = n[j] >> shift;
        HPBC_CLOCKWORK_ASSERT2(tmp[j] <= MASK);
        // normally we use (tmp & MASK), but it's redundant with tmp <= MASK
        result[j] = table[static_cast<size_t>(tmp[j])][j];
    }

    while (shift >= P) {
        for (int i=0; i<P; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for(size_t j=0; j<ARRAY_SIZE; ++j)
                result[j]= mf[j].template square<hc::LowuopsTag>(result[j]);
        }
        shift -= P;
        std::array<size_t, ARRAY_SIZE> index;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            tmp[j] = n[j] >> shift;
            index[j] = static_cast<size_t>(tmp[j]) & MASK;
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                        table[index[j]][j]);
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
    std::array<size_t, ARRAY_SIZE> index;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        index[j] = static_cast<size_t>(n[j]) & tmpmask;
        result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                        table[index[j]][j]);
    }
    return result;

  }


#ifdef __clang__
#  pragma GCC diagnostic pop
#endif
};


}} // end namespace

#endif
