// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_EXPERIMENTAL_MONTGOMERY_POW_2KARY_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_EXPERIMENTAL_MONTGOMERY_POW_2KARY_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontgomeryFormExtensions.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/count_leading_zeros.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/branchless_shift_right.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include <type_traits>
#include <cstddef>
#include <array>

namespace hurchalla { namespace experimental {



// Implementation note: this is a modified version of the 2^k-ary exponentiation
// algorithm  ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring ),
// which precalculates the even exponents as well as the normal odd exponents,
// in order to avoid two conditional branches that would exist in the main loop
// of the normal 2^k-ary algorithm.  This is particularly helpful for the array
// versions of the function further below.
//
// We use a struct with static member functions to disallow ADL.

struct experimental_montgomery_pow_2kary {


  // This function calculates pow(x, nexp), modulo the modulus of mf, and
  // returns the Montgomery domain result.
  //
  // MF must be a MontgomeryForm type (i.e. either plain MontgomeryForm, or one
  //   of its aliases like MontgomeryQuarter or MontgomeryHalf).
  // mf is an instance of MF (thus mf was constructed with a particular modulus)
  // x is the base (in Montgomery domain)
  // nexp is the exponent (it can be any unsigned integer type U).
  //
  template <class MF, typename U,
            bool USE_SLIDING_WINDOW_OPTIMIZATION = false,
            size_t TABLE_BITS = 4,
            size_t CODE_SECTION = 0,
            bool USE_SQUARING_VALUE_OPTIMIZATION = false>
  static typename MF::MontgomeryValue
  call(const MF& mf, typename MF::MontgomeryValue x, U nexp)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");
    static_assert(0 <= TABLE_BITS && TABLE_BITS < 10, "FYI you almost certainly "
        "want 2 <= TABLE_BITS <= 5.  TABLE_BITS > 0 is required.  Anything "
        "above 9 is probably a very bad idea even if it works (9+ would cause "
        "the beginning of this function to calculate 1024+ table entries!)");

    HPBC_CLOCKWORK_PRECONDITION(nexp >= 0);

    namespace hc = hurchalla;
    using V = typename MF::MontgomeryValue;
    using MFE = hc::detail::MontgomeryFormExtensions<MF, hc::LowlatencyTag>;
    using SV = typename MFE::SquaringValue;
    using std::size_t;
    U n = static_cast<U>(nexp);


if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
    // this is a masked version of pow from montgomery_pow.h
    V base = x;
    U exponent = n;

    V result;
    if (static_cast<size_t>(exponent) & 1u)
        result = base;
    else
        result = mf.getUnityValue();

    while (exponent > 1u) {
        exponent = static_cast<U>(exponent >> 1);

        base = mf.square(base);
        V tmp = mf.getUnityValue();
        // note: since we are doing masked selections, we definitely don't
        // want to use cselect_on_bit here
        tmp.template cmov<CSelectMaskedTag>(static_cast<size_t>(exponent) & 1u, base);
        result = mf.multiply(result, tmp);
    }
    return result;

} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {
        // this is a branch version of pow from montgomery_pow.h
        V base = x;
        U exponent = n;

        V result = mf.getUnityValue();
        if (static_cast<size_t>(exponent) & 1u)
            result = base;

        while (exponent > 1u) {
            exponent = static_cast<U>(exponent >> 1);
            base = mf.square(base);
            if (static_cast<size_t>(exponent) & 1u)
                result = mf.multiply(result, base);
        }
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2) {
        // this is a cmov version of pow from montgomery_pow.h
        V base = x;
        U exponent = n;

        V mont_one = mf.getUnityValue();
#ifndef HURCHALLA_MONTGOMERY_POW_2KARY_USE_CSELECT_ON_BIT
        V result = mont_one;
        result.cmov(static_cast<size_t>(exponent) & 1u, base);
#else
        V result = V::template cselect_on_bit_ne0<0>(static_cast<uint64_t>(exponent), base, mont_one);
#endif
        exponent = static_cast<U>(exponent >> 1);

        while (exponent > 0u) {
            base = mf.square(base);

#ifndef HURCHALLA_MONTGOMERY_POW_2KARY_USE_CSELECT_ON_BIT
            V tmp = mont_one;
            tmp.cmov(static_cast<size_t>(exponent) & 1u, base);
#else
            V tmp = V::template cselect_on_bit_ne0<0>(static_cast<uint64_t>(exponent), base, mont_one);
#endif
            result = mf.multiply(result, tmp);

            exponent = static_cast<U>(exponent >> 1);
        }
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 3) {
        // this is a table select adaption of CODE_SECTION 2
        V base = x;
        U exponent = n;

        V tmp[2];
        tmp[0] = mf.getUnityValue();
        tmp[1] = base;
        V result = tmp[static_cast<size_t>(exponent) & 1u];
        exponent = static_cast<U>(exponent >> 1);

        while (exponent > 0) {
            base = mf.square(base);

            tmp[1] = base;
            result = mf.multiply(result, tmp[(static_cast<size_t>(exponent) & 1u)]);

            exponent = static_cast<U>(exponent >> 1u);
        }
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 4) {
        // this is a bigger table (probably optimal size) version of CODE_SECTION 3
        V base = x;
        U exponent = n;

        V tmp[4];
        tmp[0] = mf.getUnityValue();
        tmp[1] = base;
        V result = tmp[static_cast<size_t>(exponent) & 1u];
        exponent = static_cast<U>(exponent >> 1);

        while (exponent > 0) {
            base = mf.square(base);
            V baseSqrd = mf.square(base);
            tmp[3] = mf.template multiply<LowuopsTag>(baseSqrd, base);
            tmp[1] = base;
            tmp[2] = baseSqrd;
            base = baseSqrd;
            result = mf.template multiply<LowuopsTag>(result, tmp[(static_cast<size_t>(exponent) & 3u)]);
            exponent = static_cast<U>(exponent >> 2u);
        }
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 5) {
        // standard 2k-ary table.

        constexpr int P = static_cast<int>(TABLE_BITS);
        static_assert(P > 0, "");
        constexpr size_t TABLESIZE = 1u << P;
        static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");

        V table[TABLESIZE];
        table[0] = mf.getUnityValue();   // montgomery one
        table[1] = x;
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table[2] = mf.square(x);
            table[3] = mf.multiply(table[2], x);
        }
        // if TABLESIZE is somewhat small, unroll the loop, otherwise don't unroll.
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE <= 16) {
            if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
                table[4] = mf.square(table[2]);
                table[5] = mf.multiply(table[3], table[2]);
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=6; i<TABLESIZE; i+=2) {
                    size_t j = i/2;
                    table[i] = mf.template square<LowuopsTag>(table[j]);
                    table[i+1] = mf.template multiply<LowuopsTag>(table[j + 1], table[j]);
                }
            }
        } else {
            // we should check for a ((power of 2) >= 64), but this is
            // probably adequate for our needs
            HPBC_CLOCKWORK_ASSERT(TABLESIZE % 64 == 0);
            table[4] = mf.square(table[2]);
            table[5] = mf.multiply(table[3], table[2]);
            for (size_t i=6; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table[i] = mf.template square<LowuopsTag>(table[j]);
                table[i+1] = mf.template multiply<LowuopsTag>(table[j + 1], table[j]);
            }
        }

        constexpr size_t MASK = TABLESIZE - 1;
        V result;
        if (n <= MASK) {
            result = table[static_cast<size_t>(n)];
            return result;
        }

        // count_leading_zeros returns the number of leading 0-bits in n, starting
        // at the most significant bit position. If n is 0, the result is undefined
        HPBC_CLOCKWORK_ASSERT(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        // because we returned above if (n <= MASK), we can assert the following:
        HPBC_CLOCKWORK_ASSERT(numbits > P);

        int shift = numbits - P;
        U tmp = branchless_shift_right(n, shift);
        HPBC_CLOCKWORK_ASSERT(tmp <= MASK);
        // normally we'd use (tmp & MASK), but it's redundant with tmp <= MASK
        size_t index = static_cast<size_t>(tmp);
        result = table[index];

        while (shift >= P) {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
                SV sv = MFE::getSquaringValue(mf, result);
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (shift > P && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                        sv = MFE::squareSV(mf, sv);
                        --shift;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= P);

                static_assert(P >= 1, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P-1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);
            }
            else {
                if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                    while (shift > P && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                        result = mf.square(result);
                        --shift;
                    }
                }
                HPBC_CLOCKWORK_ASSERT2(shift >= P);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<P; ++i)
                    result = mf.square(result);
            }

            shift -= P;
            index = static_cast<size_t>(branchless_shift_right(n, shift)) & MASK;
            result = mf.multiply(result, table[index]);
        }

        if (shift == 0)
            return result;
        HPBC_CLOCKWORK_ASSERT(0 < shift && shift < P);

        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        size_t tmpmask = (1u << shift) - 1;
        index = static_cast<size_t>(n) & tmpmask;
        result = mf.multiply(result, table[index]);
        return result;

} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 6) {
        // try two standard 2k-ary tables. why not?

        static_assert(TABLE_BITS > 0, "");
        constexpr size_t TABLESIZE = 1u << TABLE_BITS;
        static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");

        V table1[TABLESIZE];
        table1[0] = mf.getUnityValue();   // montgomery one
        table1[1] = x;
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table1[2] = mf.square(x);
            table1[3] = mf.multiply(table1[2], x);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            table1[4] = mf.square(table1[2]);
            table1[5] = mf.multiply(table1[3], table1[2]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=6; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table1[i] = mf.template square<LowuopsTag>(table1[j]);
                table1[i+1] = mf.template multiply<LowuopsTag>(table1[j + 1], table1[j]);
            }
        }

        V table2[TABLESIZE];
        table2[0] = mf.getUnityValue();
        table2[1] = mf.square(table1[TABLESIZE / 2]);
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table2[2] = mf.square(table2[1]);
            table2[3] = mf.multiply(table2[2], table2[1]);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            table2[4] = mf.square(table2[2]);
            table2[5] = mf.multiply(table2[3], table2[2]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=6; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table2[i] = mf.template square<LowuopsTag>(table2[j]);
                table2[i+1] = mf.template multiply<LowuopsTag>(table2[j + 1], table2[j]);
            }
        }


        constexpr size_t MASK = TABLESIZE - 1;
        constexpr int NUMBITS_MASKBIG = TABLE_BITS + TABLE_BITS;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

        V result;
        if (n <= MASKBIG) {
            size_t tmp = static_cast<size_t>(n);
            size_t loindex = tmp & MASK;
            size_t hiindex = tmp >> TABLE_BITS;
            HPBC_CLOCKWORK_ASSERT2(hiindex <= MASK);
            result = mf.multiply(table2[hiindex], table1[loindex]);
            return result;
        }

        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
        int shift = numbits - NUMBITS_MASKBIG;

        HPBC_CLOCKWORK_ASSERT2(shift > 0);
        size_t tmp = static_cast<size_t>(branchless_shift_right(n, shift));
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        size_t loindex = tmp & MASK;
        size_t hiindex = tmp >> TABLE_BITS;
        HPBC_CLOCKWORK_ASSERT2(hiindex <= MASK);
        result = mf.multiply(table2[hiindex], table1[loindex]);

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
                tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                loindex = tmp & MASK;
                hiindex = (tmp >> TABLE_BITS) & MASK;

                V val1 = mf.template multiply<LowuopsTag>(table2[hiindex], table1[loindex]);

                static_assert(NUMBITS_MASKBIG >= 1, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
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
                tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                loindex = tmp & MASK;
                hiindex = (tmp >> TABLE_BITS) & MASK;

                V val1 = mf.template multiply<LowuopsTag>(table2[hiindex], table1[loindex]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);

                result = mf.multiply(result, val1);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        tmp = static_cast<size_t>(n) & tmpmask;
        loindex = tmp & MASK;
        hiindex = tmp >> TABLE_BITS;
        HPBC_CLOCKWORK_ASSERT2(hiindex <= MASK);
        V val1 = mf.template multiply<LowuopsTag>(table2[hiindex], table1[loindex]);

        for (int i=0; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val1);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 7) {
        // three standard 2k-ary tables

        static_assert(TABLE_BITS > 0, "");
        constexpr size_t TABLESIZE = 1u << TABLE_BITS;
        static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");


        constexpr size_t MASK = TABLESIZE - 1;
        constexpr int NUMBITS_MASKBIG = 3 * TABLE_BITS;
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

        int shift = 0;
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;
        }
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);
        size_t tmp = static_cast<size_t>(branchless_shift_right(n, shift));
        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);

        size_t index1 = tmp & MASK;
        size_t index2 = (tmp >> TABLE_BITS) & MASK;
        size_t index3 = tmp >> (2*TABLE_BITS);
        HPBC_CLOCKWORK_ASSERT2(index3 <= MASK);


        V table1[TABLESIZE];
        table1[0] = mf.getUnityValue();   // montgomery one
        table1[1] = x;
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table1[2] = mf.square(x);
            table1[3] = mf.multiply(table1[2], x);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            table1[4] = mf.square(table1[2]);
            table1[5] = mf.multiply(table1[3], table1[2]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=6; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table1[i] = mf.template square<LowuopsTag>(table1[j]);
                table1[i+1] = mf.template multiply<LowuopsTag>(table1[j + 1], table1[j]);
            }
        }

        V table2[TABLESIZE];
        table2[0] = mf.getUnityValue();
        table2[1] = mf.square(table1[TABLESIZE / 2]);
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table2[2] = mf.square(table2[1]);
            table2[3] = mf.multiply(table2[2], table2[1]);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            table2[4] = mf.square(table2[2]);
            table2[5] = mf.multiply(table2[3], table2[2]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=6; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table2[i] = mf.template square<LowuopsTag>(table2[j]);
                table2[i+1] = mf.template multiply<LowuopsTag>(table2[j + 1], table2[j]);
            }
        }


        V val1 = mf.template multiply<LowuopsTag>(table2[index2], table1[index1]);


        V table3[TABLESIZE];
        table3[0] = mf.getUnityValue();
        table3[1] = mf.square(table2[TABLESIZE / 2]);
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table3[2] = mf.square(table3[1]);
            table3[3] = mf.multiply(table3[2], table3[1]);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            table3[4] = mf.square(table3[2]);
            table3[5] = mf.multiply(table3[3], table3[2]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=6; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table3[i] = mf.template square<LowuopsTag>(table3[j]);
                table3[i+1] = mf.template multiply<LowuopsTag>(table3[j + 1], table3[j]);
            }
        }


        V result = mf.multiply(table3[index3], val1);

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
                tmp = static_cast<size_t>(branchless_shift_right(n, shift));

                index1 = tmp & MASK;
                index2 = (tmp >> TABLE_BITS) & MASK;
                val1 = mf.template multiply<LowuopsTag>(table2[index2], table1[index1]);

                static_assert(NUMBITS_MASKBIG >= 3, "");
                sv = MFE::squareSV(mf, sv);
                sv = MFE::squareSV(mf, sv);

                index3 = (tmp >> (2*TABLE_BITS)) & MASK;
                V val2 = mf.template multiply<LowuopsTag>(val1, table3[index3]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=2; i<NUMBITS_MASKBIG - 1; ++i)
                    sv = MFE::squareSV(mf, sv);
                result = MFE::squareToMontgomeryValue(mf, sv);

                result = mf.multiply(result, val2);
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
                tmp = static_cast<size_t>(branchless_shift_right(n, shift));

                index1 = tmp & MASK;
                index2 = (tmp >> TABLE_BITS) & MASK;
                val1 = mf.template multiply<LowuopsTag>(table2[index2], table1[index1]);

                static_assert(NUMBITS_MASKBIG >= 2, "");
                result = mf.square(result);
                result = mf.square(result);

                index3 = (tmp >> (2*TABLE_BITS)) & MASK;
                V val2 = mf.template multiply<LowuopsTag>(val1, table3[index3]);

                HURCHALLA_REQUEST_UNROLL_LOOP for (int i=2; i<NUMBITS_MASKBIG; ++i)
                    result = mf.square(result);

                result = mf.multiply(result, val2);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        tmp = static_cast<size_t>(n) & tmpmask;

        index1 = tmp & MASK;
        index2 = (tmp >> TABLE_BITS) & MASK;
        val1 = mf.template multiply<LowuopsTag>(table2[index2], table1[index1]);

        result = mf.square(result);

        index3 = tmp >> (2*TABLE_BITS);
        HPBC_CLOCKWORK_ASSERT2(index3 <= MASK);
        V val2 = mf.template multiply<LowuopsTag>(val1, table3[index3]);

        for (int i=1; i<shift; ++i)
            result = mf.square(result);
        result = mf.multiply(result, val2);
        return result;

} else if HURCHALLA_CPP17_CONSTEXPR (8 <= CODE_SECTION && CODE_SECTION <= 20) {
        // (almost) unlimited number of tables!

        static_assert(TABLE_BITS > 0, "");
        constexpr size_t TABLESIZE = 1u << TABLE_BITS;
        static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
        constexpr size_t MASK = TABLESIZE - 1;

        constexpr size_t NUM_TABLES = CODE_SECTION - 7;
        static_assert(NUM_TABLES > 0, "");
        constexpr int NUMBITS_MASKBIG = NUM_TABLES * TABLE_BITS;
        static_assert(std::numeric_limits<size_t>::digits > NUMBITS_MASKBIG, "");
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

        std::array<std::array<V, TABLESIZE>, NUM_TABLES> table;

        table[0][0] = mf.getUnityValue();   // montgomery one
        table[0][1] = x;
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table[0][2] = mf.square(x);
            table[0][3] = mf.multiply(table[0][2], x);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=4; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table[0][i] = mf.template square<LowuopsTag>(table[0][j]);
                table[0][i+1] = mf.template multiply<LowuopsTag>(table[0][j + 1], table[0][j]);
            }
        }

        int shift = 0;
        size_t tmp = static_cast<size_t>(n);
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;

            HPBC_CLOCKWORK_ASSERT2(shift >= 0);
            tmp = static_cast<size_t>(branchless_shift_right(n, shift));
        }
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);

        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        V result = table[0][tmp & MASK];


        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k < NUM_TABLES; ++k) {
            table[k][0] = mf.getUnityValue();
            table[k][1] = mf.square(table[k - 1][TABLESIZE / 2]);
            if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
                table[k][2] = mf.square(table[k][1]);
                table[k][3] = mf.multiply(table[k][2], table[k][1]);
            }
            if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=4; i<TABLESIZE; i+=2) {
                    size_t j = i/2;
                    table[k][i] = mf.template square<LowuopsTag>(table[k][j]);
                    table[k][i+1] = mf.template multiply<LowuopsTag>(table[k][j + 1], table[k][j]);
                }
            }

            size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
            result = mf.template multiply<LowuopsTag>(table[k][index], result);

            // this part could be removed - it provides fast return when n is small.
            size_t limit_in_progress = 1u << (k * TABLE_BITS + TABLE_BITS);
            if (n < limit_in_progress)
                return result;
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
                tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                V val1 = table[0][tmp & MASK];

                static_assert(TABLE_BITS >= 1, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS - 1; ++i)
                    sv = MFE::squareSV(mf, sv);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                    tmp = tmp >> TABLE_BITS;
                    size_t index = tmp & MASK;
                    val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);

                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS; ++i)
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
                tmp = static_cast<size_t>(branchless_shift_right(n, shift));
                V val1 = table[0][tmp & MASK];

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS; ++i)
                    result = mf.square(result);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                    tmp = tmp >> TABLE_BITS;
                    size_t index = tmp & MASK;
                    val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);

                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS; ++i)
                        result = mf.square(result);
                }

                result = mf.multiply(result, val1);
            }
        }
        if (shift == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < shift && shift < NUMBITS_MASKBIG);

        size_t tmpmask = (1u << shift) - 1u;
        tmp = static_cast<size_t>(n) & tmpmask;

        V val1 = table[0][tmp & MASK];

        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
                val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);
            }

            SV sv = MFE::getSquaringValue(mf, result);
            HPBC_CLOCKWORK_ASSERT2(shift >= 1);
            for (int i=0; i<shift-1; ++i)
                sv = MFE::squareSV(mf, sv);
            result = MFE::squareToMontgomeryValue(mf, sv);
        }
        else {
            if HURCHALLA_CPP17_CONSTEXPR (NUM_TABLES == 1) {
                for (int i=0; i<shift; ++i)
                    result = mf.square(result);
            }
            else {
                static_assert(NUM_TABLES > 1, "");
                size_t index = (tmp >> TABLE_BITS) & MASK;
                val1 = mf.template multiply<LowuopsTag>(val1, table[1][index]);

                result = mf.square(result);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=2; k<NUM_TABLES; ++k) {
                    index = (tmp >> (k * TABLE_BITS)) & MASK;
                    val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);
                }

                for (int i=1; i<shift; ++i)
                    result = mf.square(result);
            }
        }

        result = mf.multiply(result, val1);
        return result;
} else if HURCHALLA_CPP17_CONSTEXPR (21 <= CODE_SECTION && CODE_SECTION <= 33) {
        // this is an optimization of the previous CODE_SECTION, to use
        // 'bits_remaining' instead of 'shift'.  It should produce more
        // efficient shifts when U is a 128 bit (or larger) integer type.

        static_assert(TABLE_BITS > 0, "");
        constexpr size_t TABLESIZE = 1u << TABLE_BITS;
        static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
        constexpr size_t MASK = TABLESIZE - 1;

        constexpr size_t NUM_TABLES = CODE_SECTION - 20;
        static_assert(NUM_TABLES > 0, "");
        constexpr int NUMBITS_MASKBIG = NUM_TABLES * TABLE_BITS;
        static_assert(std::numeric_limits<size_t>::digits > NUMBITS_MASKBIG, "");
        constexpr size_t MASKBIG = (1u << NUMBITS_MASKBIG) - 1u;

        std::array<std::array<V, TABLESIZE>, NUM_TABLES> table;

        table[0][0] = mf.getUnityValue();   // montgomery one
        table[0][1] = x;
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            table[0][2] = mf.square(x);
            table[0][3] = mf.multiply(table[0][2], x);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=4; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                table[0][i] = mf.template square<LowuopsTag>(table[0][j]);
                table[0][i+1] = mf.template multiply<LowuopsTag>(table[0][j + 1], table[0][j]);
            }
        }

        U n_orig = n;
        int shift = 0;
        size_t tmp = static_cast<size_t>(n);
        if (n > MASKBIG) {
            HPBC_CLOCKWORK_ASSERT2(n > 0);
            int leading_zeros = count_leading_zeros(n);
            int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
            HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
            shift = numbits - NUMBITS_MASKBIG;

            HPBC_CLOCKWORK_ASSERT2(shift >= 0);
            tmp = static_cast<size_t>(branchless_shift_right(n, shift));
            // this preps n ahead of time for the main loop
            n = branchless_shift_left(n, leading_zeros + NUMBITS_MASKBIG);
        }
        HPBC_CLOCKWORK_ASSERT2(shift >= 0);

        HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);
        V result = table[0][tmp & MASK];


        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k < NUM_TABLES; ++k) {
            table[k][0] = mf.getUnityValue();
            table[k][1] = mf.square(table[k - 1][TABLESIZE / 2]);
            if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
                table[k][2] = mf.square(table[k][1]);
                table[k][3] = mf.multiply(table[k][2], table[k][1]);
            }
            if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=4; i<TABLESIZE; i+=2) {
                    size_t j = i/2;
                    table[k][i] = mf.template square<LowuopsTag>(table[k][j]);
                    table[k][i+1] = mf.template multiply<LowuopsTag>(table[k][j + 1], table[k][j]);
                }
            }

            size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
            result = mf.template multiply<LowuopsTag>(table[k][index], result);

            // this part could be removed - it provides fast return when n_orig is small.
            size_t limit_in_progress = 1u << (k * TABLE_BITS + TABLE_BITS);
            if (n_orig < limit_in_progress)
                return result;
        }
        int bits_remaining = shift;


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
        // the conditional below is just to avoid a compiler warning about a
        // negative shift in the loop, even though it would never happen
        constexpr int small_shift = (digits_smaller < NUMBITS_MASKBIG)
                                     ? 0 : (digits_smaller - NUMBITS_MASKBIG);


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

                tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
                n = static_cast<U>(n << NUMBITS_MASKBIG);
                bits_remaining -= NUMBITS_MASKBIG;

                V val1 = table[0][tmp & MASK];

                static_assert(TABLE_BITS >= 1, "");
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS - 1; ++i)
                    sv = MFE::squareSV(mf, sv);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                    tmp = tmp >> TABLE_BITS;
                    size_t index = tmp & MASK;
                    val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);

                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS; ++i)
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

                tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
                n = static_cast<U>(n << NUMBITS_MASKBIG);
                bits_remaining -= NUMBITS_MASKBIG;

                V val1 = table[0][tmp & MASK];

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS; ++i)
                    result = mf.square(result);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                    tmp = tmp >> TABLE_BITS;
                    size_t index = tmp & MASK;
                    val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);

                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t i=0; i<TABLE_BITS; ++i)
                        result = mf.square(result);
                }

                result = mf.multiply(result, val1);
            }
        }
        if (bits_remaining <= 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < bits_remaining && bits_remaining < NUMBITS_MASKBIG);

//        size_t tmpmask = (1u << shift) - 1u;
//        tmp = static_cast<size_t>(n) & tmpmask;
        tmp = static_cast<size_t>(n >> high_word_shift) >> (digits_smaller - bits_remaining);

        V val1 = table[0][tmp & MASK];

        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
                val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);
            }

            SV sv = MFE::getSquaringValue(mf, result);
            HPBC_CLOCKWORK_ASSERT2(bits_remaining >= 1);
            for (int i=0; i<bits_remaining-1; ++i)
                sv = MFE::squareSV(mf, sv);
            result = MFE::squareToMontgomeryValue(mf, sv);
        }
        else {
            if HURCHALLA_CPP17_CONSTEXPR (NUM_TABLES == 1) {
                for (int i=0; i<bits_remaining; ++i)
                    result = mf.square(result);
            }
            else {
                static_assert(NUM_TABLES > 1, "");
                size_t index = (tmp >> TABLE_BITS) & MASK;
                val1 = mf.template multiply<LowuopsTag>(val1, table[1][index]);

                result = mf.square(result);

                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=2; k<NUM_TABLES; ++k) {
                    index = (tmp >> (k * TABLE_BITS)) & MASK;
                    val1 = mf.template multiply<LowuopsTag>(val1, table[k][index]);
                }

                for (int i=1; i<bits_remaining; ++i)
                    result = mf.square(result);
            }
        }

        result = mf.multiply(result, val1);
        return result;
} else {
        // this is a placeholder section til we write code to replace it
        static_assert(CODE_SECTION == 34, "");
        return mf.getUnityValue();
}
  }



#ifdef __clang__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpass-failed"
#endif

  // This is the array version of montgomery_pow_2kary.  It performs ARRAY_SIZE
  // modular exponentiations when called.  It exists as an optimization that
  // can have significantly higher throughput than the non-array version above.
  //
  template <class MF, typename U,
            size_t ARRAY_SIZE,
            size_t TABLE_BITS = 4,
            size_t CODE_SECTION = 0,
            bool USE_SQUARING_VALUE_OPTIMIZATION = false>
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf,
       const std::array<typename MF::MontgomeryValue, ARRAY_SIZE>& x,
       const std::array<U, ARRAY_SIZE>& n)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");
    static_assert(0 < TABLE_BITS && TABLE_BITS < 10, "FYI you almost certainly "
        "want 2 <= TABLE_BITS <= 5.  TABLE_BITS > 0 is required.  Anything "
        "above 9 is probably a very bad idea even if it works (9+ would cause "
        "the beginning of this function to calculate 1024+ table entries!)");

    namespace hc = hurchalla;
    using V = typename MF::MontgomeryValue;
    using MFE_LU = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;
    using SV = typename MFE_LU::SquaringValue;
    using std::size_t;

    constexpr int P = static_cast<int>(TABLE_BITS);

    // initialize the precalculation table for 2^k-ary pow algorithm
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1u << P;
    static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
    constexpr size_t MASK = TABLESIZE - 1;

    // standard 2k-ary array pow

    V table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        table[0][j] = mf[j].getUnityValue();   // montgomery one
        table[1][j] = x[j];
    }
    if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 2) {
        for (std::size_t i=2; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                table[i][j]= mf[j].template square<hc::LowuopsTag>(table[halfi][j]);
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                table[i+1][j] = mf[j].template multiply<hc::LowuopsTag>(
                                            table[halfi+1][j], table[halfi][j]);
            }
        }
    }

    U n_max = n[0];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=1; j<ARRAY_SIZE; ++j)
        n_max = (n_max < n[j]) ? n[j] : n_max;

    std::array<V, ARRAY_SIZE> result;
    if (n_max <= MASK) {
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = table[static_cast<size_t>(n[j])][j];
        return result;
    }

    // count_leading_zeros returns the number of leading 0-bits in n_max,
    // starting at the most significant bit position. If n_max is 0, the result
    // is undefined
    HPBC_CLOCKWORK_ASSERT(n_max > 0);
    int leading_zeros = count_leading_zeros(n_max);
    int numbits = ut_numeric_limits<decltype(n_max)>::digits - leading_zeros;
    // because we returned above if (n_max <= MASK), we can assert the following
    HPBC_CLOCKWORK_ASSERT(numbits > P);


    if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
        int shift = numbits - P;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            HPBC_CLOCKWORK_ASSERT(static_cast<U>(branchless_shift_right(n[j], shift)) <= MASK);
            // We don't need to 'and' with MASK, because (branchless_shift_right(n[j], shift)) <= MASK.
            size_t index = static_cast<size_t>(branchless_shift_right(n[j], shift));
            result[j] = table[index][j];
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
                size_t tmp = static_cast<size_t>(branchless_shift_right(n[j], shift));
                size_t index = tmp & MASK;
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
    } else {
        static_assert(CODE_SECTION == 1, "");
        // this is an optimization of CODE_SECTION 0, using 'bits_remaining'
        // instead of 'shift' to produce more efficient shifts when U is a
        // 128 bit (or larger) integer type.

        int bits_remaining = numbits;
        HPBC_CLOCKWORK_ASSERT2(bits_remaining > P);

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

        std::array<U, ARRAY_SIZE> n2;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            n2[j] = branchless_shift_left(n[j], leading_zeros);
            size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller - P);
            n2[j] = static_cast<U>(n2[j] << P);
            HPBC_CLOCKWORK_ASSERT2(index <= MASK);
            // normally we use (index & MASK), but it's redundant with index <= MASK
            result[j] = table[index][j];
        }
        bits_remaining -= P;

        while (bits_remaining >= P) {
            bits_remaining -= P;

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
                size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (digits_smaller-P);
                n2[j] = static_cast<U>(n2[j] << P);
                result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j],
                                                                   table[index][j]);
            }
        }
        if (bits_remaining == 0)
            return result;

        HPBC_CLOCKWORK_ASSERT2(0 < bits_remaining && bits_remaining < P);

        for (int i=0; i<bits_remaining; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf[j].template square<hc::LowuopsTag>(result[j]);
        }
        int final_shift = digits_smaller - bits_remaining;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            size_t index = static_cast<size_t>(n2[j] >> high_word_shift) >> (final_shift);
            result[j] = mf[j].template multiply<hc::LowuopsTag>(result[j], table[index][j]);
        }
        return result;
    }
  }




  // This is the "partial" array version.
  //
  // This is an alternative to calling the MontgomeryForm pow() member function
  // overload which takes an array of bases.
  // Aside from this function having a MontgomeryForm parameter (mf), it should
  // be a drop-in replacement for the MontgomeryForm pow() member function.  For
  // more information, see MontgomeryForm.h and the array version of the pow()
  // function contained within it.
  template <class MF, typename U,
            size_t ARRAY_SIZE,
            bool USE_SLIDING_WINDOW_OPTIMIZATION = false,
            size_t TABLE_BITS = 4,
            size_t CODE_SECTION = 0,
            bool USE_SQUARING_VALUE_OPTIMIZATION = false,
            class PTAG = LowuopsTag>
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const MF& mf,
       const std::array<typename MF::MontgomeryValue, ARRAY_SIZE>& x,
       U nexp)
  {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!ut_numeric_limits<U>::is_signed, "");
    static_assert(0 <= TABLE_BITS && TABLE_BITS < 10, "FYI you almost certainly "
        "want 2 <= TABLE_BITS <= 5.  TABLE_BITS > 0 is required.  Anything "
        "above 9 is probably a very bad idea even if it works (9+ would cause "
        "the beginning of this function to calculate 1024+ table entries!)");

    using V = typename MF::MontgomeryValue;
    using MFE_LU = hurchalla::detail::MontgomeryFormExtensions<MF, PTAG>;
    using SV = typename MFE_LU::SquaringValue;
    using std::size_t;
    U n = nexp;

    constexpr int P = static_cast<int>(TABLE_BITS);

    // initialize the precalculation table for 2^k-ary pow algorithm
    static_assert(P > 0, "");
    constexpr size_t TABLESIZE = 1u << P;
    static_assert(TABLESIZE >= 2 && TABLESIZE % 2 == 0, "");
    constexpr size_t MASK = TABLESIZE - 1;


if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 0) {
    // this is adapted from arraypow_cond_branch_unrolled() in montgomery_pow.h

    std::array<V, ARRAY_SIZE> bases = x;
    U exponent = n;

    std::array<V, ARRAY_SIZE> result;
    if (static_cast<size_t>(exponent) & 1u) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = bases[j];
    } else {
        V mont_one = mf.getUnityValue();
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mont_one;
    }

    while (exponent > 1u) {
        exponent = static_cast<U>(exponent >> 1);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            bases[j] = mf.template square<PTAG>(bases[j]);
        if (static_cast<size_t>(exponent) & 1u) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf.template multiply<PTAG>(result[j], bases[j]);
        }
    }
    return result;

} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 1) {

    V table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        table[0][j] = mf.getUnityValue();   // montgomery one
        table[1][j] = x[j];
    }
    if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
        for (std::size_t i=2; i<TABLESIZE; i+=2) {
            std::size_t halfi = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i][j] = mf.template square<PTAG>(table[halfi][j]);
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                table[i+1][j] = mf.template multiply<PTAG>(table[halfi+1][j], table[halfi][j]);
        }
    }

    std::array<V, ARRAY_SIZE> result;
    if (n <= MASK) {
        size_t index = static_cast<size_t>(n);
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = table[index][j];
        return result;
    }

    // count_leading_zeros returns the number of leading 0-bits in n, starting
    // at the most significant bit position. If n is 0, the result is undefined
    HPBC_CLOCKWORK_ASSERT(n > 0);
    int leading_zeros = count_leading_zeros(n);
    int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
    // because we returned above if (n <= MASK), we can assert the following:
    HPBC_CLOCKWORK_ASSERT(numbits > P);

    int shift = numbits - P;
    size_t index = static_cast<size_t>(branchless_shift_right(n, shift));
    HPBC_CLOCKWORK_ASSERT(index <= MASK);
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        // normally we'd use (index & MASK), but it's redundant with index <= MASK
        result[j] = table[index][j];
    }


    while (shift >= P) {
        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            std::array<SV, ARRAY_SIZE> sv;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                sv[j] = MFE_LU::getSquaringValue(mf, result[j]);

            if (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > P && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf, sv[j]);
                    --shift;
                }
            }

            static_assert(P > 0, "");
            for (int i=0; i<P - 1; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::squareSV(mf, sv[j]);
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = MFE_LU::squareToMontgomeryValue(mf, sv[j]);
        }
        else {
            if (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > P && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf.template square<PTAG>(result[j]);
                    --shift;
                }
            }

            for (int i=0; i<P; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = mf.template square<PTAG>(result[j]);
            }
        }

        shift -= P;
        index = static_cast<size_t>(branchless_shift_right(n, shift)) & MASK;
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            result[j] = mf.template multiply<PTAG>(result[j], table[index][j]);
        }
    }

    if (shift == 0)
        return result;
    HPBC_CLOCKWORK_ASSERT(0 < shift && shift < P);

    for (int i=0; i<shift; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mf.template square<PTAG>(result[j]);
    }
    size_t tmpmask = (1u << shift) - 1;
    index = static_cast<size_t>(n) & tmpmask;
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        result[j] = mf.template multiply<PTAG>(result[j], table[index][j]);
    }
    return result;
} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 2) {

    // This CODE_SECTION optimizes table initialization to skip the high even
    // elements of the table.  The while loop does clever cmovs to avoid ever
    // needing to use those (uninitialized) high even elements, without using
    // branches or any extra multiplies.
    //
    // This method could also be used to create another version of the full
    // array 2k-ary pow function, but since this CODE_SECTION hasn't benchmarked
    // any better than CODE_SECTION 1, and since this method is somewhat less
    // suited to full array 2k-ary pow than it is to partial array pow, it seems
    // very unlikely that this method would improve the full array pow.

    static_assert(TABLESIZE >= 8, "");

    constexpr int NUMBITS_MASK_SMALL = P - 1;
    static_assert(NUMBITS_MASK_SMALL >= 0, "");
    constexpr size_t MASK_SMALL = (1u << NUMBITS_MASK_SMALL) - 1u;
    static_assert(MASK_SMALL == (TABLESIZE/2) - 1, "");

    V table[TABLESIZE][ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
        table[0][j] = mf.getUnityValue();   // montgomery one
        table[1][j] = x[j];
    }
    if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 2) {
        constexpr size_t HALFSIZE = TABLESIZE/2;
        for (std::size_t i=2; i<HALFSIZE; i+=2) {
            std::size_t halfi = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                table[i][j]= mf.template square<PTAG>(table[halfi][j]);
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                table[i+1][j] = mf.template multiply<PTAG>(
                                            table[halfi+1][j], table[halfi][j]);
            }
        }
        constexpr size_t QUARTERSIZE = TABLESIZE/4;
        for (std::size_t i=1; i<HALFSIZE; i+=2) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                table[HALFSIZE + i][j] = mf.template multiply<PTAG>(
                              table[QUARTERSIZE + i][j], table[QUARTERSIZE][j]);
            }
        }
    }

    std::array<V, ARRAY_SIZE> result;
    if (n <= MASK_SMALL) {
        HPBC_CLOCKWORK_ASSERT(n <= (TABLESIZE/2) - 1);
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = table[static_cast<size_t>(n)][j];
        return result;
    }

    // count_leading_zeros returns the number of leading 0-bits in n,
    // starting at the most significant bit position. If n is 0, the result
    // is undefined
    HPBC_CLOCKWORK_ASSERT(n > 0);
    int leading_zeros = count_leading_zeros(n);
    int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
    // because we returned above if (n <= MASK_SMALL), we can assert the following
    HPBC_CLOCKWORK_ASSERT(numbits >= P);

    int shift = numbits - P;
    HPBC_CLOCKWORK_ASSERT(shift >= 0);
    HPBC_CLOCKWORK_ASSERT((branchless_shift_right(n, shift)) <= MASK);
    // due to above assert, we don't need to 'and' with MASK
    size_t index = static_cast<size_t>(branchless_shift_right(n, shift));

    // because the highest set bit of n is by definition a 1, we know
    HPBC_CLOCKWORK_ASSERT((index >> (P-1)) == 1u);  // and thus
    HPBC_CLOCKWORK_ASSERT(index >= TABLESIZE/2);

    size_t index1 =  2 * index - (TABLESIZE - 1);
    size_t index2 = (TABLESIZE - 1) - index;
    HPBC_CLOCKWORK_ASSERT(index1 % 2 == 1);
    HPBC_CLOCKWORK_ASSERT(index2 < TABLESIZE/2);
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
        result[j] = mf.template multiply<PTAG>(table[index1][j], table[index2][j]);


    while (shift >= P) {
        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            std::array<SV, ARRAY_SIZE> sv;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                sv[j] = MFE_LU::getSquaringValue(mf, result[j]);

            if (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > P && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        sv[j] = MFE_LU::squareSV(mf, sv[j]);
                    --shift;
                }
            }

            // we do P-1 squarings instead of P squarings
            static_assert(P > 1, "");
            for (int i=0; i<P - 2; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    sv[j] = MFE_LU::squareSV(mf, sv[j]);
            }
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = MFE_LU::squareToMontgomeryValue(mf, sv[j]);
        }
        else {
            if (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (shift > P && (static_cast<size_t>(branchless_shift_right(n, shift-1)) & 1u) == 0) {
                    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                        result[j] = mf.template square<PTAG>(result[j]);
                    --shift;
                }
            }

            // we do P-1 squarings instead of P squarings
            static_assert(P > 0, "");
            for (int i=0; i<P - 1; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                    result[j] = mf.template square<PTAG>(result[j]);
            }
        }

        shift -= P;
        index = static_cast<size_t>(branchless_shift_right(n, shift)) & MASK;

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
#ifndef HURCHALLA_MONTGOMERY_POW_2KARY_USE_CSELECT_ON_BIT
            V tmp = result[j];
            tmp.cmov((index % 2 == 0), table[index/2][j]);
#else
            V tmp = V::template cselect_on_bit_eq0<0>(static_cast<uint64_t>(index), table[index/2][j], result[j]);
#endif
            result[j] = mf.template multiply<PTAG>(tmp, result[j]);
        }

        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
#ifndef HURCHALLA_MONTGOMERY_POW_2KARY_USE_CSELECT_ON_BIT
            V tmp = table[index][j];
            tmp.cmov((index % 2 == 0), result[j]);
#else
            V tmp = V::template cselect_on_bit_eq0<0>(static_cast<uint64_t>(index), result[j], table[index][j]);
#endif
            result[j] = mf.template multiply<PTAG>(tmp, result[j]);
        }
    }

    if (shift == 0)
        return result;
    HPBC_CLOCKWORK_ASSERT(0 < shift && shift < P);

    for (int i=0; i<shift; ++i) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mf.template square<PTAG>(result[j]);
    }

    size_t tmpmask = (1u << shift) - 1;
    HPBC_CLOCKWORK_ASSERT(tmpmask <= MASK_SMALL);
    index = static_cast<size_t>(n) & tmpmask;
    HPBC_CLOCKWORK_ASSERT(index < TABLESIZE/2);
    HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
        result[j] = mf.template multiply<PTAG>(result[j],table[index][j]);

    return result;

} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 3) {
    // this is adapted from arraypow_cmov() in montgomery_pow.h

    std::array<V, ARRAY_SIZE> bases = x;
    U exponent = n;

    std::array<V, ARRAY_SIZE> result;
    if (static_cast<size_t>(exponent) & 1u) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = bases[j];
    } else {
        V mont_one = mf.getUnityValue();
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mont_one;
    }

    while (exponent > 1u) {
        exponent = static_cast<U>(exponent >> 1);
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            bases[j] = mf.template square<PTAG>(bases[j]);

        V mont_one = mf.getUnityValue();
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
# ifndef HURCHALLA_MONTGOMERY_POW_2KARY_USE_CSELECT_ON_BIT
            V tmp = mont_one;
            tmp.cmov(static_cast<size_t>(exponent) & 1u, bases[j]);
# else
            V tmp = V::template cselect_on_bit_ne0<0>(
                           static_cast<uint64_t>(exponent), bases[j], mont_one);
# endif
            result[j] = mf.template multiply<PTAG>(result[j], tmp);
        }
    }
    return result;

} else if HURCHALLA_CPP17_CONSTEXPR (CODE_SECTION == 4) {
    // this is adapted from arraypow_masked() in montgomery_pow.h

    std::array<V, ARRAY_SIZE> bases = x;
    U exponent = n;

    std::array<V, ARRAY_SIZE> result;
    if (static_cast<size_t>(exponent) & 1u) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = bases[j];
    } else {
        V mont_one = mf.getUnityValue();
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mont_one;
    }

    while (exponent > 1u) {
        exponent = static_cast<U>(exponent >> 1);

        V mont_one = mf.getUnityValue();
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            bases[j] = mf.template square<PTAG>(bases[j]);
            V tmp = mont_one;
            // note: since we are doing masked selections, we definitely don't
            // want to use cselect_on_bit here
            tmp.template cmov<CSelectMaskedTag>(static_cast<size_t>(exponent) & 1u, bases[j]);
            result[j] = mf.template multiply<PTAG>(result[j], tmp);
        }
    }
    return result;

} else {
    static_assert(CODE_SECTION == 5, "");
    // this is adapted from arraypow_cond_branch() in montgomery_pow.h

    std::array<V, ARRAY_SIZE> bases = x;
    U exponent = n;

    std::array<V, ARRAY_SIZE> result;
    if (static_cast<size_t>(exponent) & 1u) {
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = bases[j];
    } else {
        V mont_one = mf.getUnityValue();
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            result[j] = mont_one;
    }

    while (exponent > 1u) {
        exponent = static_cast<U>(exponent >> 1);
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            bases[j] = mf.template square<PTAG>(bases[j]);
        if (static_cast<size_t>(exponent) & 1u) {
            for (size_t j=0; j<ARRAY_SIZE; ++j)
                result[j] = mf.template multiply<PTAG>(result[j], bases[j]);
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
