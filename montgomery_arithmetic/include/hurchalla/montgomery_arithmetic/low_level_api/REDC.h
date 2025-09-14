// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_LLAPI_REDC_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_LLAPI_REDC_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/ImplRedc.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <type_traits>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla {


// This file is the API for the alternate REDC algorithm described at
// https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/README_REDC.md
// This alternate version of the REDC algorithm differs in small but important
// ways from Peter Montgomery's original 1985 paper "Modular multiplication
// without trial division".  From the point of view of a caller, the most
// important distinction is that this version requires the positive inverse for
// one of its arguments rather than the negative inverse (which was required by
// the original/traditional REDC algorithm).  We provide the alternate version
// instead of the traditional version, because it improves efficiency both in
// terms of latency and number of instructions.  See README_REDC.md for the
// details.
//
// For discussion purposes below, let the unlimited precision constant R
// represent  R = 1<<(ut::ut_numeric_limits<T>::digits).  For example, if T is
// uint64_t, then R = 1<<64.


//   REDC_standard() returns the standard and normally expected value from REDC,
// which is the least residue (modulo the modulus).  In other words,
// REDC_standard() guarantees that 0 <= return_value < modulus.
//   When calling REDC_standard() you must specify a PTAG type - use one of the
// structs from  modular_arithmetic/detail/optimization_tag_structs.h  and see
// optimization_tag_structs.h for the benefits of different PTAGs.
template <typename T, class PTAG> HURCHALLA_FORCE_INLINE
T REDC_standard(T u_hi, T u_lo, T n, T inv_n, PTAG)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");

    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);
    using P = typename safely_promote_unsigned<T>::type;
    // verify that  n * inv_n ≡ 1 (mod R)
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);  // verify that (u_hi*R + u_lo) < n*R

    T result = detail::RedcStandard<T>::call(u_hi, u_lo, n, inv_n, PTAG());

    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
}


// REDC_incomplete() is "incomplete" in that this function does not perform the
// final subtraction and does not conditionally add the modulus to that
// difference, both of which would be needed to obtain a completed REDC result.
// Instead, it provides the minuend and subtrahend, allowing the caller to
// perform the eventual final subtraction and (usually) the conditional add.
// The caller can perform these final steps in whatever way is most suitable to
// its needs.
// The reason you might wish to use this function is that it can provide better
// performance than the standard REDC in some situations.
// ---
// As an example of an optimization that REDC_incomplete() allows, we can use it
// to optimize Montgomery multiplication when the modulus 'n' is less than R/4.
// For such a case, we can unconditionally add 'n' to the difference of minuend
// and subtrahend (REDC_incomplete provides these values), and we can then use
// this sum "as-is" as an input to a Montgomery multiplication.  There is no
// need to use an extra instruction (or more) to make this a conditional add of
// the modulus that would minimize this input value to the least residue (i.e.
// to 0 <= sum < n), because Montgomery multiplication with inputs x and y only
// requires that  u = x*y < n*R,  and we will always be able to satisfy this
// requirement when using an unconditionally added sum, provided that the
// modulus 'n' is less than R/4.  For details see section 5 from
// "Montgomery's Multiplication Technique: How to Make It Smaller and Faster"
// https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps
// You can see also MontyQuarterRange.h in this library, which is a class
// that requires n < R/4 and is optimized in this way.
// Another example of a different optimization enabled by REDC_incomplete is
// MontyHalfrange.h.


// When calling REDC_incomplete() you must specify a PTAG type - use one of the
// structs from  modular_arithmetic/detail/optimization_tag_structs.h  and see
// optimization_tag_structs.h for the benefits of different PTAGs.
template <typename T, class PTAG> HURCHALLA_FORCE_INLINE
void REDC_incomplete(T& minuend, T& subtrahend, T u_hi, T u_lo, T n, T inv_n, PTAG)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");

    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);
    using P = typename safely_promote_unsigned<T>::type;
    // verify that  n * inv_n ≡ 1 (mod R)
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);  // verify that (u_hi*R + u_lo) < n*R

    detail::RedcIncomplete::call(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());

    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T diff = static_cast<T>(minuend - subtrahend);
        T finalized_result = (minuend < subtrahend) ? static_cast<T>(diff + n) : diff;
        HPBC_CLOCKWORK_POSTCONDITION2(finalized_result ==
                      ::hurchalla::REDC_standard(u_hi, u_lo, n, inv_n, PTAG()));
        HPBC_CLOCKWORK_POSTCONDITION2(0 <= finalized_result && finalized_result < n);
    }
    // If  n < R/2,  then  0 < diff + n < 2*n
    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T diff = static_cast<T>(minuend - subtrahend);
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_CLOCKWORK_POSTCONDITION2((n < Rdiv2) ? (0 < static_cast<T>(diff + n)) &&
                                (static_cast<T>(diff + n) < 2*n) : true);
    }
}


#if 0
// Obsolete version, but it should work fine if enabled by setting the #if to 1.
//
//   When calling REDC_incomplete() you must specify a PTAG type - use one of
// the structs from  modular_arithmetic/detail/optimization_tag_structs.h  and
// see optimization_tag_structs.h for the benefits of different PTAGs.
//
template <typename T, class PTAG> HURCHALLA_FORCE_INLINE
T REDC_incomplete(bool& isNegative, T u_hi, T u_lo, T n, T inv_n, PTAG)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");

    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);
    using P = typename safely_promote_unsigned<T>::type;
    // verify that  n * inv_n ≡ 1 (mod R)
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);  // verify that (u_hi*R + u_lo) < n*R

    T minuend, subtrahend;
    detail::RedcIncomplete::call(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    T result = static_cast<T>(minuend - subtrahend);
    isNegative = minuend < subtrahend;

    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T finalized_result = (isNegative) ? static_cast<T>(result + n) : result;
        HPBC_CLOCKWORK_POSTCONDITION2(finalized_result ==
                      ::hurchalla::REDC_standard(u_hi, u_lo, n, inv_n, PTAG()));
        HPBC_CLOCKWORK_POSTCONDITION2(0 <= finalized_result && finalized_result < n);
    }
    // If  n < R/2,  then  0 < result + n < 2*n
    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_CLOCKWORK_POSTCONDITION2((n < Rdiv2) ? (0 < static_cast<T>(result + n)) &&
                                (static_cast<T>(result + n) < 2*n) : true);
    }
    return result;
}
#endif



// This version of REDC_incomplete calculates the final difference (of minuend
// and subtrahend) and returns it, but it does not add the modulus to the
// difference in any way before returning.  Since this return value might be
// positive or negative (with no indication which), this is an incomplete REDC.
// This function can be useful when either you know ahead of time that the
// result will be negative or non-negative, or when it doesn't matter whether or
// not the result will be negative.
//
// When calling REDC_incomplete() you must specify a PTAG type - use one of the
// structs from  modular_arithmetic/detail/optimization_tag_structs.h  and see
// optimization_tag_structs.h for the benefits of different PTAGs.
template <typename T, class PTAG> HURCHALLA_FORCE_INLINE
T REDC_incomplete(T u_hi, T u_lo, T n, T inv_n, PTAG)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");

    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);
    using P = typename safely_promote_unsigned<T>::type;
    // verify that  n * inv_n ≡ 1 (mod R)
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);  // verify that (u_hi*R + u_lo) < n*R

    T result = detail::RedcIncomplete::call(u_hi, u_lo, n, inv_n, PTAG());

    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        namespace hc = ::hurchalla;
        T complete_result = hc::REDC_standard(u_hi, u_lo, n, inv_n, PTAG());
        HPBC_CLOCKWORK_POSTCONDITION2(result == complete_result ||
                                 static_cast<T>(result + n) == complete_result);
    }
    // If  n < R/2,  then  0 < result + n < 2*n
    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_CLOCKWORK_POSTCONDITION2((n < Rdiv2) ? (0 < static_cast<T>(result + n)) &&
                                (static_cast<T>(result + n) < 2*n) : true);
    }
    return result;
}


} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
