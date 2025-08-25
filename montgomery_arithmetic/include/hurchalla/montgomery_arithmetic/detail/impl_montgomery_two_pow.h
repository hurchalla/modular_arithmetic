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
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>
#include <cstddef>
#include <array>
#include <utility>

namespace hurchalla { namespace detail {


// Implementation note: this is a highly modified version of the 2^k-ary
// algorithm  ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring ),
// taking advantage of the fact that the base is always 2.  It also
// precalculates the even exponents as well as the normal odd exponents, in
// order to avoid the conditional branches that would exist in the main loop of
// the normal 2^k-ary algorithm.  This is particularly helpful for the array
// version of this function further below.
//
// Minor note: we use a struct with static member functions to disallow ADL.
struct impl_montgomery_two_pow {

  template <typename U>
  static constexpr int floor_log2(U x)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");
    // x > 0 required, but C++11 constexpr function won't allow this check
    //HPBC_PRECONDITION2(x > 0);
    return (x <= 1) ? 0 : 1 + impl_montgomery_two_pow::floor_log2(x >> 1);
  }


  // Calculate pow(2, n), modulo the modulus of mf, and return the result in
  // montgomeryform representation.
  template <bool USE_SLIDING_WINDOW_OPTIMIZATION,
            size_t TABLE_BITS,
            size_t CODE_SECTION,
            class MF, typename U>
  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    namespace hc = hurchalla;
    using V = typename MF::MontgomeryValue;
    using MFE = hc::detail::MontgomeryFormExtensions<MF, hc::LowlatencyTag>;
    using std::size_t;

    using RU = typename MFE::RU;
    constexpr int digitsRU = hurchalla::ut_numeric_limits<RU>::digits;
    constexpr int P2 = floor_log2(static_cast<unsigned int>(digitsRU));
    constexpr size_t MASK = (1u << P2) - 1u;

    // note: The unused template parameter TABLE_BITS exists simply to maintain
    // correspondence with experimental_montgomery_two_pow.cpp.

    // note: This code is a copy of experimental_montgomery_two_pow.cpp's
    // non-array call(), for TABLE_BITS == 0, CODE_SECTION 1, 2, and 3 (you can
    // find code section 0 there too, purposely not included here).
    HPBC_ASSERT2(TABLE_BITS == 0 && (CODE_SECTION == 1 || CODE_SECTION == 2 ||
                                     CODE_SECTION == 3));

    if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
        if (n <= MASK) {
            size_t loindex = static_cast<size_t>(n);
            RU num = static_cast<RU>(static_cast<RU>(1) << loindex);
            return MFE::convertInExtended(mf, num);
        }
        RU magicValue = MFE::getMagicValue(mf);

        HPBC_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_ASSERT2(numbits >= (P2 + 1));

        int shift = numbits - (P2 + 1);
        HPBC_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_ASSERT2(tmp <= 2u*MASK + 1u);
        // Bit P2 of tmp was the leading bit, so it should always be set.
        HPBC_ASSERT2(((tmp >> P2) & 1u) == 1u);
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
                HPBC_ASSERT2(shift >= (P2 + 1));

                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                num = static_cast<RU>(static_cast<RU>(1) << loindex);
                V val1 = MFE::convertInExtended_aTimesR(mf, num, magicValue);
                HPBC_ASSERT2(((tmp >> P2) & 1u) == 1u);
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

        HPBC_ASSERT2(0 < shift && shift < (P2 + 1));

        size_t tmpmask = (1u << shift) - 1u;
        size_t index = static_cast<size_t>(n) & tmpmask;
        RU num2 = static_cast<RU>(static_cast<RU>(1) << index);
        V tableVal = MFE::convertInExtended(mf, num2);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
    }
    else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2) {
        if (n <= MASK)
            return MFE::twoPowLimited(mf, static_cast<size_t>(n));
        HPBC_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_ASSERT2(numbits > P2);

        int shift = numbits - P2;
        size_t index = static_cast<size_t>(n >> shift);
        HPBC_ASSERT2(index <= MASK);
        V result = MFE::twoPowLimited(mf, index);
        while (shift >= P2) {
            if (USE_SLIDING_WINDOW_OPTIMIZATION) {
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
        HPBC_ASSERT2(0 < shift && shift < P2);

        size_t tmpmask = (1u << shift) - 1u;
        index = static_cast<size_t>(n) & tmpmask;
        V tableVal = MFE::twoPowLimited(mf, index);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
    } else {        // CODE_SECTION 3
        if (n <= MASK) {
            size_t loindex = static_cast<size_t>(n);
            return MFE::twoPowLimited(mf, loindex);
        }
        RU magicValue = MFE::getMagicValue(mf);

        HPBC_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_ASSERT2(numbits >= (P2 + 1));

        int shift = numbits - (P2 + 1);
        HPBC_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(n >> shift);
        HPBC_ASSERT2(tmp <= 2u*MASK + 1u);
        // Bit P2 of tmp was the leading bit, so it should always be set.
        HPBC_ASSERT2(((tmp >> P2) & 1u) == 1u);
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
                HPBC_ASSERT2(shift >= (P2 + 1));

                shift -= (P2 + 1);
                tmp = static_cast<size_t>(n >> shift);
                loindex = tmp & MASK;
                V val1 = MFE::RTimesTwoPowLimited(mf, loindex, magicValue);
                HPBC_ASSERT2(((tmp >> P2) & 1u) == 1u);
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

        HPBC_ASSERT2(0 < shift && shift < (P2 + 1));

        size_t tmpmask = (1u << shift) - 1u;
        size_t index = static_cast<size_t>(n) & tmpmask;
        V tableVal = MFE::twoPowLimited(mf, index);
        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, tableVal);
        return result;
    }
  }





#ifdef __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpass-failed"
#endif

  // Array version of montgomery two pow
  template <size_t TABLE_BITS,
            size_t CODE_SECTION,
            class MF, typename U, size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");

    namespace hc = hurchalla;
    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;
    using MFE = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;
    using std::size_t;

    namespace hc = hurchalla;

    U n_max = n[0];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=1; j<ARRAY_SIZE; ++j) {
        if (n_max < n[j])
            n_max = n[j];
    }

    // note: The unused template parameter TABLE_BITS exists simply to maintain
    // correspondence with experimental_montgomery_two_pow.cpp.

    // note: This code is a copy of experimental_montgomery_two_pow.cpp's
    // array call(), for TABLE_BITS == 0, CODE_SECTION 0 and 2 (you can find
    // code section 1 there too, which is purposely not included here).

    HPBC_ASSERT2(TABLE_BITS == 0 && (CODE_SECTION == 0 || CODE_SECTION == 2));

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

        HPBC_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_ASSERT2(numbits > P2);

        int shift = numbits - P2;
        std::array<V, ARRAY_SIZE> result;
        std::array<size_t, ARRAY_SIZE> tmp;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            tmp[j] = static_cast<size_t>(n[j] >> shift);
            HPBC_ASSERT2(tmp[j] <= MASK);
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
        HPBC_ASSERT2(0 < shift && shift < P2);
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
    else {      // CODE_SECTION == 2
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
        HPBC_ASSERT2(n_max > 0);
        int leading_zeros = count_leading_zeros(n_max);
        int numbits= ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
        HPBC_ASSERT2(numbits > 3);
        int shift = numbits - 3;

        std::array<V, ARRAY_SIZE> result;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n[j] >> shift);
            HPBC_ASSERT2(index <= MASK);
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


#ifdef __clang__
#  pragma GCC diagnostic pop
#endif
};


}} // end namespace

#endif
