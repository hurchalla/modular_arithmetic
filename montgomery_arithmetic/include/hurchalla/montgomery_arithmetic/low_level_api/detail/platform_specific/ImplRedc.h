// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_REDC_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_REDC_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace detail {


// Montgomery REDC algorithm
//
// This file implements the REDC algorithm as described at
// https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/README_REDC.md
// This is an alternate version of the REDC algorithm, that differs in small but
// important ways from Peter Montgomery's original 1985 paper "Modular
// multiplication without trial division".  From the point of view of a caller,
// the most important distinction is that this version requires the positive
// inverse for one of its arguments rather than the negative inverse (which was
// required by the original/traditional REDC algorithm).  For our purposes, the
// most important distinction is that this alternate version is a more efficient
// algorithm both for latency and number of instructions.  See README_REDC.md
// for the details.
// Note that the description in README_REDC.md uses a variable name "T", despite
// the fact that in C++ "T" is conventionally reserved for use as a template
// parameter name.   This is done for consistency with most all descriptions of
// Montgomery multiplication/REDC, including Montgomery's 1985 paper, the
// Wikipedia article https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
// and many more.

// In all functions below, we will use the variable name "u" (and u_hi and u_lo)
// in place of the algorithm description's variable name "T", since as already
// stated, in C++ "T" is de-facto a reserved name for a template parameter.  We
// will use "T" in its normal C++ sense, as a template parameter name - this has
// no relation to "T" in the algorithm description.
// We also use "n" instead of "N", and "inv_n" instead of "N^(-1)" (N with the
// superscript -1).  The constant "R" remains the same, and represents the value
// R = 1<<(ut::ut_numeric_limits<T>::digits).  As an example, if T is uint64_t,
// then R = 1<<64.


// The name "RedcIncomplete" reflects the fact that the call() function's return
// value is not finalized to the least residue mod n (i.e. the return value does
// not necessarily satisfy 0 <= return_value < n).  See the starred* comment
// under Postcondition #1 for more info.
// Minor note: uses static member function to disallow ADL.
struct RedcIncomplete {
  template <typename T>
  static HURCHALLA_FORCE_INLINE T call(bool& ovf, T u_hi, T u_lo, T n, T inv_n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");

    // For casts, we want to use types that are protected from surprises and
    // undefined behavior due to the unsigned integral promotion rules in C++.
    // https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/
    using P = typename safely_promote_unsigned<T>::type;
    static_assert(ut_numeric_limits<P>::is_integer, "");
    static_assert(!(ut_numeric_limits<P>::is_signed), "");
    static_assert(ut_numeric_limits<P>::is_modulo, "");

    // Precondition #1:  We require the precondition  u < n*R.  Or elaborated,
    // u == u_hi*R + u_lo < n*R.
    // If u_hi < n:  then u_hi+1 <= n, and u_hi*R + R <= n*R.  Since u_lo < R,
    //   u == u_hi*R + u_lo < u_hi*R + R <= n*R.  We would have  u < n*R, and
    //   so u_hi < n is always sufficient to satisfy the precondition.
    // If u_hi >= n:  then u_hi*R >= n*R, and u == u_hi*R + u_lo >= n*R, which
    //   fails the precondition.
    // Thus u_hi < n is sufficient and necessary to satisfy the precondition.
    HPBC_PRECONDITION2(u_hi < n);

    // assert(n * inv_n ≡ 1 (mod R))
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    // compute  m = (u * inv_n) % R
    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));

    T mn_lo;
    T mn_hi = ::hurchalla::unsigned_multiply_to_hilo_product(mn_lo, m, n);

    // mn = m*n.  Since m = (u_lo*inv_n)%R, we know m < R, and thus  mn < R*n.
    // Therefore mn == mn_hi*R + mn_lo < R*n, and mn_hi*R < R*n - mn_lo <= R*n,
    // and thus  mn_hi < n.
        // *** Assertion #1 ***
    HPBC_ASSERT2(mn_hi < n);

    // The REDC algorithm from README_REDC.md assures us that (u - mn) is
    // divisible by R.  Compute (u - mn)/R  (and we can note that a negative
    // result in C++ will wrap around to a large value, very similar to two's
    // complement representation of a negative value) :
    T t_hi = static_cast<T>(u_hi - mn_hi);   // t_hi = (u_hi - mn_hi) mod R

    ovf = (u_hi < mn_hi);    // tells us if the subtraction wrapped/overflowed

    // We do not need to explicitly perform the low part subtraction
    // (u_lo - mn_lo), because the REDC algorithm guarantees
    // (u_lo - mn_lo) mod R == 0.  Since 0 <= u_lo < R and 0 <= mn_lo < R, this
    // means that u_lo == mn_lo, and thus (u_lo - mn_lo) will never generate a
    // borrow/carry.  We will simply ignore this low part subtraction.
        // *** Assertion #2 ***
    HPBC_ASSERT2(u_lo == mn_lo);

    // Since u_hi and u_lo are type T (which is unsigned) variables, both
    // u_hi >= 0 and u_lo >= 0, and thus  u = u_hi*R + u_lo >= 0.  Along with
    // the precondition u < n*R, we therefore know that  0 <= u < n*R.
    // Since m is type T we know m >= 0, and since n is type T we know n >= 0.
    // Thus mn >= 0.  Subtracting mn from all parts of 0 <= u < n*R, we have
    // -mn <= u - mn < n*R - mn.  Since we just showed mn >= 0, we know
    // -mn <= 0 and thus  n*R - mn <= n*R.  Assertion #1 states that mn < n*R,
    // and thus we know  -n*R < -mn.  Therefore
    // -n*R < -mn <= u - mn < n*R - mn <= n*R.  Simplified, -n*R < u - mn < n*R.
    // Since (u - mn) is divisible by R (see README_REDC.md for proof), and
    // R > 0, we have  -n < (u - mn)/R < n.
    // Since u = u_hi*R + u_lo, and mn = mn_hi*R + mn_lo,  we know
    // u - mn == u_hi*R + u_lo - mn_hi*R - mn_lo.  Assertion #2 states that
    // u_lo == mn_lo, and so we have
    // u - mn == u_hi*R - mn_hi*R == (u_hi - mn_hi)*R,  and thus
    // (u - mn)/R == u_hi - mn_hi.  Therefore,  -n < u_hi - mn_hi < n.
    //
    // All this translates into
    // Postcondition #1
    // ----------------
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T finalized_result = (ovf) ? static_cast<T>(t_hi + n) : t_hi;
        HPBC_POSTCONDITION2(finalized_result < n);
    }
    // * Aside from this postcondition, we do not actually compute the finalized
    // least residual mod n result, because some Montgomery Forms are
    // constrained in ways that allow a simpler and more efficient computation
    // of the finalized result.  For example, in some forms the input u_hi (and
    // the return value) is allowed to occupy the range 0 <= u_hi < 2*n, which
    // lets us change the conditional add of n at the end of the REDC into an
    // unconditional add of n.

    // Postcondition #2:  If  n < R/2,  then  0 < t_hi + n < 2*n
    // ---------------------------------------------------------
    // We already showed  -n < u_hi - mn_hi < n.  Adding n to all parts we get
    // 0 < u_hi - mn_hi + n < 2*n.  Although this is true regardless of the size
    // of n, we can only test this postcondition when n < R/2 (any larger value
    // of n would overflow on 2*n).
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (0 < static_cast<T>(t_hi + n)) &&
                                       (static_cast<T>(t_hi + n) < 2*n) : true);
    }

    // return the non-finalized result
    return t_hi;
  }
};



template <typename T>
struct DefaultRedcStandard
{
  static_assert(ut_numeric_limits<T>::is_integer, "");
  static_assert(!(ut_numeric_limits<T>::is_signed), "");
  static_assert(ut_numeric_limits<T>::is_modulo, "");

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n)
  {
    bool ovf;
    T result = RedcIncomplete::call(ovf, u_hi, u_lo, n, inv_n);
    // By RedcIncomplete::call()'s Postcondition #1, we get
    // T final_result = (ovf) ? static_cast<T>(result + n) : result;
#if defined(HURCHALLA_AVOID_CSELECT)
    T mask = static_cast<T>(-static_cast<T>(ovf));
    T final_result = static_cast<T>(result + (mask & n));
#else
      // final_result = (ovf) ? static_cast<T>(result+n) : result
    T final_result = ::hurchalla::conditional_select(
                                         ovf, static_cast<T>(result+n), result);
#endif
    HPBC_POSTCONDITION2(final_result < n);
    return final_result;
  }
};


// primary template
template <typename T>
struct RedcStandard
{
  static_assert(ut_numeric_limits<T>::is_integer, "");
  static_assert(!(ut_numeric_limits<T>::is_signed), "");
  static_assert(ut_numeric_limits<T>::is_modulo, "");

  // For PTAGs see optimization_tag_structs.h
  template <class PTAG>
  static HURCHALLA_FORCE_INLINE T call(T u_hi, T u_lo, T n, T inv_n, PTAG)
  {
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n);
    HPBC_POSTCONDITION2(result < n);
    return result;
  }
};


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// specialization for uint64_t (for x86_64)
template <>
struct RedcStandard<std::uint64_t>
{
  using T = std::uint64_t;

  // This function should have: cycles latency 9, fused uops 7
  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowlatencyTag)
  {
    using P = typename safely_promote_unsigned<T>::type;
    // This implementation is based closely on
    // DefaultRedcStandard<uint64_t>::call and the RedcIncomplete::call that it
    // in turn calls.
    // Thus the algorithm should be correct for the same reasons given there.
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // RedcIncomplete's call, u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));
    T mn_lo;
    T mn_hi = ::hurchalla::unsigned_multiply_to_hilo_product(mn_lo, m, n);
    HPBC_ASSERT2(mn_hi < n);
    T reg = u_hi + n;
    T uhi = u_hi;
    __asm__ (
        "subq %[mnhi], %[reg] \n\t"     /* reg = u_hi + n - mn_hi */
        "subq %[mnhi], %[uhi] \n\t"     /* t_hi = u_hi - mn_hi */
        "cmovaeq %[uhi], %[reg] \n\t"   /* reg = (u_hi >= mn_hi) ? t_hi : reg */
        : [reg]"+&r"(reg), [uhi]"+&r"(uhi)
        : [mnhi]"r"(mn_hi)
        : "cc");
    T result = reg;
    HPBC_ASSERT2(result == DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n));
    HPBC_POSTCONDITION2(result < n);
    return result;
  }

  // This function should have: cycles latency 10, fused uops 6
  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowuopsTag)
  {
    // Calling DefaultRedcStandard::call will give us optimal code (we're
    // relying upon modular_subtract_prereduced_inputs() being optimized for low
    // uops - which it is, at least at the time of writing this)
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n);
    HPBC_POSTCONDITION2(result < n);
    return result;
  }
};
#endif   // (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) ||
         //  defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) &&
         // defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
