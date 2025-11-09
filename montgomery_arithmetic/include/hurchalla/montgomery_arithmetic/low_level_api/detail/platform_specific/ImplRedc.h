// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_REDC_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_REDC_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/unsigned_multiply_to_hi_product.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
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
// R = (UP)1 << ut::ut_numeric_limits<T>::digits, where UP is a conceptually
// unlimited precision unsigned integer type.  As an example, if T is uint64_t,
// then R = (UP)1 << 64.


struct RedcIncomplete {
// The name "RedcIncomplete" reflects the fact that this function does not
// perform the final subtraction needed to obtain a completed REDC result.
// Instead, it provides the minuend and subtrahend, allowing the caller to
// perform the eventual final subtraction.  The caller can perform the final
// subtraction in whatever way most suitable to its needs.
// Minor note: uses static member functions to disallow ADL.

  // This function calculates the minuend and subtrahend for the complete REDC,
  // such that the completed REDC =
  //        (minuend < subtrahend) ? minuend-subtrahend + n : minuend-subtrahend
  template <typename T, class PTAG>
  static HURCHALLA_FORCE_INLINE
  void call(T& minuend, T& subtrahend, T u_hi, T u_lo, T n, T inv_n, PTAG)
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
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);

    // assert(n * inv_n ≡ 1 (mod R))
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    // compute  m = (u * inv_n) % R
    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));

    //T mn_lo;
    // since we now skip calculating mn_lo (it's not used anyway), by calling
    // unsigned_multiply_to_hi_product(), mn_lo can be considered to be an
    // implied variable for explanation purposes, rather than a programming
    // variable.
    //T mn_hi = ::hurchalla::unsigned_multiply_to_hilo_product(mn_lo, m, n);
    T mn_hi = ::hurchalla::unsigned_multiply_to_hi_product(m, n);

    // mn = m*n.  Since m = (u_lo*inv_n)%R, we know m < R, and thus  mn < R*n.
    // Therefore mn == mn_hi*R + mn_lo < R*n, and mn_hi*R < R*n - mn_lo <= R*n,
    // and thus  mn_hi < n.
        // *** Assertion #1 ***
    HPBC_CLOCKWORK_ASSERT2(mn_hi < n);

    // The REDC algorithm from README_REDC.md assures us that (u - mn) is
    // divisible by R.  Compute (u - mn)/R   (note that a negative result in C++
    // will wrap around, usually resulting in a large value, using the two's
    // complement representation of a negative value) :
    //T t_hi = static_cast<T>(u_hi - mn_hi);   // t_hi = (u_hi - mn_hi) mod R

    //ovf = (u_hi < mn_hi);    // tells us if the subtraction wrapped/overflowed

    // However, in this function we are interested in only setting the minuend
    // and subtrahend.  We leave it to the caller to perform the eventual
    // subtraction.  Since (u - mn) is divisible by R, we disregard the low
    // words of the minuend and subtrahend, and thus we set the minuend to u_hi,
    // set the subtrahend to mn_hi.

    minuend = u_hi;
    subtrahend = mn_hi;

    // We do not need to explicitly provide the low parts for subtraction
    // (u_lo - mn_lo), because the REDC algorithm guarantees
    // (u_lo - mn_lo) mod R == 0.  Since 0 <= u_lo < R and 0 <= mn_lo < R, this
    // means that u_lo == mn_lo, and thus (u_lo - mn_lo) can never generate a
    // borrow/carry.
    // We simply disregard the low words of the minuend and subtrahend, since
    // they are equal to each other.
        // *** Assertion #2 ***
    // since we skip calculting the unused variable mn_lo, we can't actually
    // call this assert anymore, even though it would be true.
    //HPBC_CLOCKWORK_ASSERT2(u_lo == mn_lo);

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
    // Since  minuend == u_hi  and  subtrahend == mn_hi,  we have
    // -n < minuend - subtrahend < n.
    //
    // All this translates into
    // Postcondition #1
    // ----------------
    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T difference = static_cast<T>(minuend - subtrahend);
        bool ovf = (minuend < subtrahend);
        T finalized_result = ovf ? static_cast<T>(difference + n) : difference;
        HPBC_CLOCKWORK_POSTCONDITION2(finalized_result < n);
    }
    // * Aside from this postcondition, we do not actually compute the finalized
    // least residual mod n result, because some Montgomery Forms are
    // constrained in ways that allow a simpler and more efficient computation
    // of the finalized result.  For example, in some forms the input u_hi (and
    // the return value) is allowed to occupy the range 0 <= u_hi < 2*n, which
    // lets us change the conditional-add of n at the end of a complete REDC
    // into an unconditional add of n.

    // Postcondition #2:  If  n < R/2,  then  0 < minuend - subtrahend + n < 2*n
    // ---------------------------------------------------------
    // We already showed  -n < minuend - subtrahend < n.  Adding n to all parts
    // we get  0 < minuend - subtrahend + n < 2*n.  Although this is true
    // regardless of the size of n, we can only test this postcondition when
    // n < R/2 (any larger value of n would overflow on 2*n).
    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T t_hi = static_cast<T>(minuend - subtrahend);
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_CLOCKWORK_POSTCONDITION2((n < Rdiv2) ? (0 < static_cast<T>(t_hi + n)) &&
                                       (static_cast<T>(t_hi + n) < 2*n) : true);
    }

    // We have as precondition u_hi < n, and we set minuend == u_hi, thus
    // minuend < n.  We know from assertion #1 that  mn_hi < n, and we set
    // subtrahend = mn_hi, and thus  subtrahend < n.  So we know
    //
    // Postcondition #3:  minuend < n  and  subtrahend < n.
    HPBC_CLOCKWORK_POSTCONDITION2(minuend < n && subtrahend < n);
  }


  // This version of the function returns the (semi) final REDC result without
  // the adjustment of adding the modulus when that result is negative.  This
  // provides maximum efficiency when it doesn't matter whether the final REDC
  // result is positive or negative.
  template <typename T, class PTAG>
  static HURCHALLA_FORCE_INLINE T call(T u_hi, T u_lo, T n, T inv_n, PTAG)
  {
    T minuend;
    T subtrahend;
    call<T, PTAG>(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    HPBC_CLOCKWORK_ASSERT2(minuend < n && subtrahend < n);
    T t_hi = static_cast<T>(minuend - subtrahend);
    return t_hi;
  }



#if (HURCHALLA_COMPILER_HAS_UINT128_T())
  // The performance for these __uint128_t versions on m2 is excellent, so long
  // as throughput is needed rather than low latency.
  // Performance benefits on x64 are similar to ARM64 (m2) - these are much
  // faster than the ordinary versions when using LowuopsTag (for throughput),
  // and slower when using LowlatencyTag.


  // Calculates the minuend and subtrahend of the REDC, such that the finalized
  // REDC = (minuend < subtrahend) ? minuend-subtrahend + n : minuend-subtrahend
  //
  // This algorithm is very loosely based on the multiprecision REDC at
  // https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#Montgomery_arithmetic_on_multiprecision_integers
  static HURCHALLA_FORCE_INLINE
  void call(__uint128_t& minuend, __uint128_t& subtrahend, __uint128_t u_hi,
            __uint128_t u_lo, __uint128_t n, __uint128_t inv_n, LowuopsTag)
  {
    using T = __uint128_t;
    using TH = uint64_t;    // THalf bits

    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);
    HPBC_CLOCKWORK_PRECONDITION2(n * inv_n == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    constexpr int HALF_BITS = ut_numeric_limits<TH>::digits;
    TH n0 = static_cast<TH>(n);
    TH n1 = static_cast<TH>(n >> HALF_BITS);

    TH u0 = static_cast<TH>(u_lo);
    TH u1 = static_cast<TH>(u_lo >> HALF_BITS);

    TH invn0 = static_cast<TH>(inv_n);


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) && \
    defined(HURCHALLA_TARGET_ISA_ARM_64) && !defined(_MSC_VER)

# if 0
// this #if section corresponds to the first #if section of C++ code
// lower in this function.
// On M2, this section is about 10% slower than the #else (for both clang and
// gcc), and so it is disabled.  In theory it has an advantage, by having one
// fewer mult than the #else section, but this seems to be overcome by its
// disadvantages of needing four or five additional instructions (not needed by
// the #else section) to make a potential negative value positive during
// calculations (done by the instructions that conditionally set moz and then
// add moz).

    TH u2 = static_cast<TH>(u_hi);
    TH u3 = static_cast<TH>(u_hi >> HALF_BITS);
    TH m, tmp;
    __asm__ ("mul %[m], %[u0], %[invn0] \n\t"
             "umulh %[u0], %[m], %[n0] \n\t"     /* u0 = mnA_1 */
             "mul %[tmp], %[m], %[n1] \n\t"      /* tmp = mnA_1_part2 */
             "umulh %[m], %[m], %[n1] \n\t"      /* m = mnA_2 */
             "adds %[u0], %[tmp], %[u0] \n\t"    /* mnA_1 += mnA_1_part2 */
             "cinc %[m], %[m], hs \n\t"          /* mnA_2 += carry */
             "subs %[u1], %[u1], %[u0] \n\t"     /* u1 = v1 = u1 - mnA_1 */
             "mul %[u0], %[u1], %[invn0] \n\t"   /* u0 = mB = mullo(v1, invn0) */

             "sbcs %[u1], %[u2], %[m] \n\t"      /* u1 = v2 = u2 - mnA_2 - borrow */
             "sbcs %[tmp], %[u3], xzr \n\t"      /* tmp = v3 = u3 - borrow */
#   if 1
             /* on an M2 this perfs better than the #else by ~2% */
             /* note: on M2, csel and csetm run on exec units 1-3; 'and' runs on units 1-6 */
             "csetm %[invn0], lo \n\t"           /* invn0 = mask = -borrow */
             "and %[u2], %[n0], %[invn0] \n\t"   /* u2 = moz0 = n0 & mask */
             "and %[u3], %[n1], %[invn0] \n\t"   /* u3 = moz1 = n1 & mask */
#   else
             "csel %[u2], %[n0], xzr, lo \n\t"   /* u2 = moz0 = (u3<borrow) ? n0 : 0 */
             "csel %[u3], %[n1], xzr, lo \n\t"   /* u3 = moz1 = (u3<borrow) ? n1 : 0 */
#   endif
             "adds %[u1], %[u1], %[u2] \n\t"     /* v2 += moz0 */
             "adcs %[tmp], %[tmp], %[u3] \n\t"   /* v3 += moz1 */

             "umulh %[m], %[u0], %[n0] \n\t"     /* m = mnB_2 = mulhi(mB, n0) */
             "mul %[u2], %[u0], %[n1] \n\t"      /* u2 = mnB_2_part2 = mullo(mB, n1) */
             "umulh %[u0], %[u0], %[n1] \n\t"    /* u0 = mnB_3 = mulhi(mB, n1) */
             "adds %[m], %[u2], %[m] \n\t"       /* mnB_2 += mnB_2_part2 */
             "cinc %[u0], %[u0], hs \n\t"        /* mnB_3 += carry */
             : [m]"=&r"(m), [u0]"+&r"(u0), [invn0]"+&r"(invn0),
               [n0]"+&r"(n0), [tmp]"=&r"(tmp), [n1]"+&r"(n1),
               [u1]"+&r"(u1), [u2]"+&r"(u2), [u3]"+&r"(u3)
             :
             : "cc");
    minuend = (static_cast<T>(tmp) << HALF_BITS) | u1;
    subtrahend = (static_cast<T>(u0) << HALF_BITS) | m;
# else
    TH m, tmp;
    __asm__ ("mul %[m], %[u0], %[invn0] \n\t"
             "umulh %[u0], %[m], %[n0] \n\t"     /* u0 = mnA_1 */
             "mul %[tmp], %[m], %[n1] \n\t"      /* tmp = mnA_1_part2 */
             "umulh %[m], %[m], %[n1] \n\t"      /* m = mnA_2 */
             "adds %[u0], %[tmp], %[u0] \n\t"    /* mnA_1 += mnA_1_part2 */
             "cinc %[m], %[m], hs \n\t"          /* mnA_2 += carry */
             "sub %[u1], %[u1], %[u0] \n\t"      /* u1 = v1 = u1 - mnA_1 */
             "mul %[tmp], %[u1], %[invn0] \n\t"  /* tmp = mB = mullo(v1, invn0) */

             "mul %[u1], %[tmp], %[n0] \n\t"     /* u1 = mnB_1 = mullo(mB, n0) */
             "umulh %[n0], %[tmp], %[n0] \n\t"   /* n0 = mnB_2 = mulhi(mB, n0) */
             "mul %[invn0], %[tmp], %[n1] \n\t"  /* invn0 = mnB_2_part2 = mullo(mB, n1) */
             "umulh %[tmp], %[tmp], %[n1] \n\t"  /* tmp = mnB_3 = mulhi(mB, n1) */
             "adds %[n0], %[invn0], %[n0] \n\t"  /* mnB_2 += mnB_2_part2 */
             "cinc %[tmp], %[tmp], hs \n\t"      /* mnB_3 += carry */

             "adds %[u0], %[u0], %[u1] \n\t"     /* u0 = dummy = mnA_1 + mnB_1 */
             "adcs %[u1], %[n0], %[m] \n\t"      /* u1 = sum2 = mnB_2 + mnA_2 + carry */
             "cinc %[tmp], %[tmp], hs \n\t"      /* tmp = sum3 = mnB_3 += carry */
             : [m]"=&r"(m), [u0]"+&r"(u0), [invn0]"+&r"(invn0),
               [n0]"+&r"(n0), [tmp]"=&r"(tmp), [n1]"+&r"(n1),
               [u1]"+&r"(u1)
             :
             : "cc");
    minuend = u_hi;
    subtrahend = (static_cast<T>(tmp) << HALF_BITS) | u1;
# endif

#elif (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

    TH tmp = u0;
    TH rrax = n0;
    TH rrdx;
    __asm__ ("imulq %[invn0], %[tmp] \n\t" /* tmp = mA = u0 * inv_n */
             "mulq %[tmp] \n\t"            /* rdx:rax = mnA_10 = rax * mA (rax == n0); high-order bits of the product in rdx */
             "movq %[tmp], %%rax \n\t"     /* rax = mA */
             "movq %%rdx, %[tmp] \n\t"     /* tmp = mnA_1 */
             "mulq %[n1] \n\t"             /* rdx:rax = mnA_21 = n1 * mA */
             "addq %%rax, %[tmp] \n\t"     /* tmp = mnA_1 += mnA_1_part2 */

             "movq %[n0], %%rax \n\t"      /* rax = n0_original */
             "movq %%rdx, %[n0] \n\t"      /* n0 = mnA_2 */

             "adcq $0, %[n0] \n\t"         /* mnA_2 += carry */
             "subq %[tmp], %[u1] \n\t"     /* u1 = v1 = u1 - mnA_1 */
             "imulq %[u1], %[invn0] \n\t"  /* invn0 = mB = v1 * invn0 */

             "mulq %[invn0] \n\t"          /* rdx:rax = mnB_21 = n0_original * mB */
             "movq %%rax, %[u1] \n\t"      /* u1 = mnB_1 */

             "movq %[invn0], %%rax \n\t"   /* rax = mB */
             "movq %%rdx, %[invn0] \n\t"   /* invn0 = mnB_2 */

             "mulq %[n1] \n\t"             /* rdx:rax = mnB_32 = n1 * mB */
             "addq %%rax, %[invn0] \n\t"   /* invn0 = mnB_2 += mnB_2_part2 */
             "adcq $0, %%rdx \n\t"         /* rdx = mnB_3 += carry */

             "addq %[u1], %[tmp] \n\t"     /* tmp = dummy = mnA_1 + mnB_1 */
             "adcq %[invn0], %[n0] \n\t"   /* n0 = sum2 = mnB_2 + mnA_2 + carry */
             "adcq $0, %%rdx \n\t"         /* rdx = sum3 = mnB_3 += carry */
             : [invn0]"+&r"(invn0), "+&a"(rrax), "=&d"(rrdx),
               [tmp]"+&r"(tmp), [u1]"+&r"(u1), [n0]"+&r"(n0)
             : [n1]"r"(n1)
             : "cc");
    minuend = u_hi;
    subtrahend = (static_cast<T>(rrdx) << HALF_BITS) | n0;

#else  // not using inline-asm


    TH mA = u0 * invn0;

    T mnA_10 = static_cast<T>(mA) * n0;      // mnA_10 <= (R-1)*(R-1) == R^2 - 2R + 1
    T mnA_21 = static_cast<T>(mA) * n1;      // mnA_21 <= R^2 - 2R + 1
    mnA_21 = mnA_21 + (mnA_10 >> HALF_BITS); // mnA_21 <= R^2 - 2R + 1 + floor((R^2 - 2R + 1)/R)
                                             // mnA_21 <= R^2 - 2R + 1 + R - 2 == R^2 - R - 1
    // sanity check:  (R^2 - 1) * (R - 1) == R^3 - R^2 - R + 1
    //                floor((R^3 - R^2 - R + 1)/R) == R^2 - R - 1.  Looks good; this is same as what we got for max possible mnA_21.
    HPBC_CLOCKWORK_ASSERT2(u0 == static_cast<TH>(mnA_10));

    TH mnA_1 = static_cast<TH>(mnA_21);
    TH v1 = u1 - mnA_1;

    TH mB = v1 * invn0;

# if 0
    T u_21 = (u_hi << HALF_BITS) | u1;
    // we skip the lowest part of the subtract (u_lo_lo - mn_lo_lo), because
    // it's implicitly a zero result and generates no borrow.

    T v_21 = u_21 - mnA_21;
    TH u3 = static_cast<TH>(u_hi >> HALF_BITS);
    // the following can go negative/underflow, which is why we set up moz
    TH v3 = u3 - (u_21 < mnA_21);
        // T moz = (u3 < (u_21 < mnA_21)) ? n : 0;
    T moz = ::hurchalla::conditional_select((u3 < (u_21 < mnA_21)), n, static_cast<__uint128_t>(0));

    T v_32 = (static_cast<T>(v3) << HALF_BITS) | (v_21 >> HALF_BITS);
    v_32 = v_32 + moz;

    T mnB_21 = static_cast<T>(mB) * n0;      // mnB_21 <= (R-1)*(R-1) == R^2 - 2R + 1
    T mnB_32 = static_cast<T>(mB) * n1;      // mnB_32 <= R^2 - 2R + 1
    mnB_32 = mnB_32 + (mnB_21 >> HALF_BITS); // mnB_32 <= R^2 - R    (see mnA_21 calculation)
    HPBC_CLOCKWORK_ASSERT2(static_cast<TH>(v_21) == static_cast<TH>(mnB_21));

    // in our function outputs of minuend and subtrahend, we don't include
    // v1 or mnB_1, because  v1 - mnB_1  is implicitly a zero result and generate
    // generates no borrow.

//    T t_hi = v_32 - mnB_32;
//    ovf = (v_32 < mnB_32);
    minuend = v_32;
    subtrahend = mnB_32;
# else
    T mnB_21 = static_cast<T>(mB) * n0;     // mnB_21 <= (R-1)*(R-1) == R^2 - 2R + 1
    T mnB_32 = static_cast<T>(mB) * n1;      // mnB_32 <= R^2 - 2R + 1
    mnB_32 = mnB_32 + (mnB_21 >> HALF_BITS);  // mnB_32 <= R^2 - R - 1    (see mnA_21 calculation)

    TH mnB_1 = static_cast<TH>(mnB_21);
    HPBC_CLOCKWORK_ASSERT2(v1 == mnB_1);

    // t = ((u - mnA)/R - mnB)/R
    // t = ((u_upper3_words - mnA_upper3_words) - mnB)/R
    // t = (u_upper3_words - (mnA_upper3_words + mnB))/R

    // we know that  u_upper3_words - mnA_upper3_words ≡ mnB  (mod R)
    // and thus      u_upper3_words ≡ mnA_upper3_words + mnB  (mod R)

    TH dummy = mnA_1 + mnB_1;
    HPBC_CLOCKWORK_ASSERT2(u1 == dummy);
    // Thus we don't need to perform  u_lo_lo - dummy  (or more precisely, we
    // don't need to provide u_lo_lo and dummy as part of the function's outputs
    // minuend and subtrahend),  because the result is implicitly zero and the
    // sub doesn't generate a borrow.

    // We do need to make sure we account for the potential carry from the
    // addition that created dummy (the carry is why dummy exists).

    TH mnA_2 = static_cast<TH>(mnA_21 >> HALF_BITS); // mnA_2 <= floor((R^2 - R - 1)/R) == R - 2
                                                     // mnA_2 <= R - 2
    T sum_hi = mnB_32 + mnA_2 + (dummy < mnB_1);     // sum_hi <= R^2 - R - 1 + R - 2 + 1 == R^2 - 2
    // unless I've made some bad mistakes, the calculation of sum_hi can not overflow.

//    T t_hi = u_hi - sum_hi;
//    ovf = (u_hi < sum_hi);
    minuend = u_hi;
    subtrahend = sum_hi;
# endif

#endif  // choice of inline-asm vs not inline-asm

    HPBC_CLOCKWORK_POSTCONDITION2(minuend < n && subtrahend < n);

    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T t_hi = minuend - subtrahend;
        T finalized_result = (minuend < subtrahend) ? t_hi + n : t_hi;
        HPBC_CLOCKWORK_POSTCONDITION2(finalized_result < n);
        T minu2, subt2;
        call<T, LowuopsTag>(minu2, subt2, u_hi, u_lo, n, inv_n, LowuopsTag());
        T answer = minu2 - subt2;
        T finalized_answer = (minu2 < subt2) ? answer + n : answer;
        HPBC_CLOCKWORK_POSTCONDITION2(finalized_result == finalized_answer);
    }
    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T t_hi = minuend - subtrahend;
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_CLOCKWORK_POSTCONDITION2((n < Rdiv2) ? (0 < (t_hi + n)) &&
                                       ((t_hi + n) < 2*n) : true);
    }
  }



  // we can implement the above algorithm more straightforwardly and more
  // efficiently here, since we return the final subtraction result while
  // making no distinction between a positive or negative result.
  static HURCHALLA_FORCE_INLINE
  __uint128_t call(__uint128_t u_hi, __uint128_t u_lo, __uint128_t n, __uint128_t inv_n, LowuopsTag)
  {
    using T = __uint128_t;
    using TH = uint64_t;    // THalf bits

    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);
    HPBC_CLOCKWORK_PRECONDITION2(n * inv_n == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    constexpr int HALF_BITS = ut_numeric_limits<TH>::digits;
    TH n0 = static_cast<TH>(n);
    TH n1 = static_cast<TH>(n >> HALF_BITS);
    TH invn0 = static_cast<TH>(inv_n);

    TH u0 = static_cast<TH>(u_lo);
    TH u1 = static_cast<TH>(u_lo >> HALF_BITS);


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) && \
    defined(HURCHALLA_TARGET_ISA_ARM_64) && !defined(_MSC_VER)

    TH u2 = static_cast<TH>(u_hi);
    TH u3 = static_cast<TH>(u_hi >> HALF_BITS);

    TH m, tmp;
    __asm__ ("mul %[m], %[u0], %[invn0] \n\t"
             "umulh %[u0], %[m], %[n0] \n\t"     /* u0 = mnA_1 */
             "mul %[tmp], %[m], %[n1] \n\t"
             "umulh %[m], %[m], %[n1] \n\t"      /* m = mnA_2 */
             "adds %[u0], %[tmp], %[u0] \n\t"    /* mnA_1 += tmp */
             "cinc %[m], %[m], hs \n\t"          /* mnA_2 += carry */
             "subs %[u1], %[u1], %[u0] \n\t"     /* u1 = v1 = u1 - mnA_1 */
             "mul %[u0], %[u1], %[invn0] \n\t"   /* u0 = mB = mullo(v1, invn0) */
             "sbcs %[u1], %[u2], %[m] \n\t"      /* u1 = v2 = u2 - mnA_2 - borrow */
             "sbc %[tmp], %[u3], xzr \n\t"       /* tmp = v3 = u3 - borrow */
             "umulh %[m], %[u0], %[n0] \n\t"     /* m = mnB_2 = mulhi(mB, n0) */
             "mul %[u2], %[u0], %[n1] \n\t"      /* u2 = mnB_2_part2 =mullo(mB, n1) */
             "umulh %[u0], %[u0], %[n1] \n\t"    /* u0 = mnB_3 = mulhi(mB, n1) */
             "adds %[m], %[u2], %[m] \n\t"       /* mnB_2 += mnB_2_part2 */
             "cinc %[u0], %[u0], hs \n\t"        /* mnB_3 += carry */
             "subs %[u1], %[u1], %[m] \n\t"      /* t2 = v2 - mnB_2 */
             "sbc %[tmp], %[tmp], %[u0] \n\t"    /* t3 = v3 - mnB_3 - borrow */
#if 1
/* strangely I've observed both gcc and clang timing about 1% faster here when
   I declare inputs as in/out params, rather than declaring them input params
   (which is much more normal).  It's difficult to say if that's likely to be
   consistent, or it was a fluke with a particular benchmark.  I've also
   observed gcc go badly wrong on benchmark timings with this asm here, where
   using this asm made the benchmark take 1.35x longer than not using the asm
   (I suspect a bug in gcc- maybe it was redundantly using this asm code
   twice?).  When I used the #else clause (where inputs are declared normally)
   the perf improved to only a 1.1x slowdown, which shows the #else can be
   better at times, even though the perf was still quite bad.
   In theory... we would expect that declaring inputs as in/out params when
   they are truly input-only (i.e. unchanged in the asm) would not improve
   perf... */
             : [m]"=&r"(m), [u0]"+&r"(u0), [invn0]"+&r"(invn0),
               [n0]"+&r"(n0), [tmp]"=&r"(tmp), [n1]"+&r"(n1),
               [u1]"+&r"(u1), [u2]"+&r"(u2), [u3]"+&r"(u3)
             :
#else
             : [m]"=&r"(m), [u0]"+&r"(u0), [tmp]"=&r"(tmp),
               [u1]"+&r"(u1), [u2]"+&r"(u2)
             : [invn0]"r"(invn0), [n0]"r"(n0), [n1]"r"(n1), [u3]"r"(u3)
#endif
             : "cc");

    T t_hi = (static_cast<T>(tmp) << HALF_BITS) | u1;

#elif (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

    TH u2 = static_cast<TH>(u_hi);
    TH u3 = static_cast<TH>(u_hi >> HALF_BITS);

    TH tmp = u0;
    TH rrax = n0;
    TH rrdx;
    __asm__ ("imulq %[invn0], %[tmp] \n\t" /* tmp = mA = u0 * inv_n */
             "mulq %[tmp] \n\t"            /* rdx:rax = mnA_10 = rax * mA (rax == n0); high-order bits of the product in rdx */
             "movq %[tmp], %%rax \n\t"     /* rax = mA */
             "movq %%rdx, %[tmp] \n\t"     /* tmp = mnA_1 */
             "mulq %[n1] \n\t"             /* rdx:rax = mnA_21 = n1 * mA */
             "addq %%rax, %[tmp] \n\t"     /* tmp = mnA_1 += mnA_1_part2 */

             "movq %[n0], %%rax \n\t"      /* rax = n0_original */

             "adcq $0, %%rdx \n\t"         /* mnA_2 += carry */
             "subq %[tmp], %[u1] \n\t"     /* u1 = v1 = u1 - mnA_1 */
             "sbbq %%rdx, %[u2] \n\t"      /* u2 = v2 = u2 - mnA_2 - borrow */
             "sbbq $0, %[u3] \n\t"         /* u3 = v3 = u3 - borrow */

             "imulq %[invn0], %[u1] \n\t"  /* u1 = mB = v1 * invn0 */

             "mulq %[u1] \n\t"             /* rdx:rax = mnB_21 = n0_original * mB */
             "movq %[u1], %%rax \n\t"      /* rax = mB */
             "movq %%rdx, %[u1] \n\t"      /* invn0 = mnB_2 */
             "mulq %[n1] \n\t"             /* rdx:rax = mnB_32 = n1 * mB */
             "addq %%rax, %[u1] \n\t"      /* invn0 = mnB_2 += mnB_2_part2 */
             "adcq $0, %%rdx \n\t"         /* rdx = mnB_3 += carry */

             "subq %[u1], %[u2] \n\t"      /* t2 = v2 - mnB_2 */
             "sbbq %%rdx, %[u3] \n\t"      /* t3 = v3 - mnB_3 - borrow */
             : [invn0]"+&r"(invn0), "+&a"(rrax), "=&d"(rrdx), [tmp]"+&r"(tmp),
               [u1]"+&r"(u1), [u2]"+&r"(u2), [u3]"+&r"(u3)
             : [n0]"r"(n0), [n1]"r"(n1)
             : "cc");
    T t_hi = (static_cast<T>(u3) << HALF_BITS) | u2;

#else
// no inline-asm

    TH mA = u0 * invn0;

    T mnA_10 = static_cast<T>(mA) * n0;
    T mnA_21 = static_cast<T>(mA) * n1;
    mnA_21 = mnA_21 + (mnA_10 >> HALF_BITS);
    HPBC_CLOCKWORK_ASSERT2(u0 == static_cast<TH>(mnA_10));

    TH mnA_1 = static_cast<TH>(mnA_21);
    TH v1 = u1 - mnA_1;
    TH mB = v1 * invn0;

# if defined(__clang__)
    // clang produces better assembly with this section (gcc does very badly on it!)
    T u_21 = (u_hi << HALF_BITS) | u1;
    T v_21 = u_21 - mnA_21;
    TH u3 = static_cast<TH>(u_hi >> HALF_BITS);
    // note: this subtraction can produce a negative number (interpreting the
    // result as a signed integer) or unsigned underflow (intepreting the
    // result as unsigned).  Both are fine for this function.
    TH v3 = u3 - (u_21 < mnA_21);
    T v_32 = (static_cast<T>(v3) << HALF_BITS) | (v_21 >> HALF_BITS);
# else
    TH mnA_2 = static_cast<TH>(mnA_21 >> HALF_BITS);
    // note: this subtraction can produce a negative number (interpreting the
    // result as a signed integer) or unsigned underflow (intepreting the
    // result as unsigned).  Both are fine for this function.
    T v_32 = (u_hi - static_cast<T>(mnA_2)) - (u1 < mnA_1);
# endif

    T mnB_21 = static_cast<T>(mB) * n0;
    T mnB_32 = static_cast<T>(mB) * n1;
    mnB_32 = mnB_32 + (mnB_21 >> HALF_BITS);
    HPBC_CLOCKWORK_ASSERT2(v1 == static_cast<TH>(mnB_21));

    // this subtraction can go negative/underflow too, which is no problem for
    // this function
    T t_hi = v_32 - mnB_32;

#endif


    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minuend, subt;
        call<T, LowuopsTag>(minuend, subt, u_hi, u_lo, n, inv_n, LowuopsTag());
        T answer = minuend - subt;
        answer = (minuend < subt) ? answer + n : answer;
        HPBC_CLOCKWORK_POSTCONDITION2(t_hi == answer || t_hi + n == answer);
    }
    if (HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_CLOCKWORK_POSTCONDITION2((n < Rdiv2) ? (0 < (t_hi + n)) &&
                                       ((t_hi + n) < 2*n) : true);
    }
    return t_hi;
  }

#endif   // endif uint128 supported



#if 0
  // Thus function is obsolete, but is expected to work fine if you enable it.
  //
  // This is an obsolete version of RedcIncomplete:call, which returns a
  // slightly more complete REDC result - it performs the final subtraction,
  // though it does not add in the modulus when the difference is negative.
  // Instead it sets  ovf = true  when the difference is negative, and sets
  // ovf = false  when the difference is positive, allowing the caller to make
  // the final adjustments.
  //
  template <typename T, class PTAG>
  static HURCHALLA_FORCE_INLINE T call(bool& ovf, T u_hi, T u_lo, T n, T inv_n, PTAG)
  {
    T minuend, subtrahend;
    call(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    HPBC_CLOCKWORK_ASSERT2(minuend < n && subtrahend < n);
    T t_hi = static_cast<T>(minuend - subtrahend);
    ovf = (minuend < subtrahend);
    return t_hi;
  }
#endif

};






#if 1
// This is the old code, which is well proven.  The #else section should in
// theory be preferable, but it needs to be tested for correctness and speed.


template <typename T>
struct DefaultRedcStandard
{
  static_assert(ut_numeric_limits<T>::is_integer, "");
  static_assert(!(ut_numeric_limits<T>::is_signed), "");
  static_assert(ut_numeric_limits<T>::is_modulo, "");

  // For PTAGs see optimization_tag_structs.h
  template <class PTAG>
  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, PTAG)
  {
    namespace hc = ::hurchalla;

    T minuend, subtrahend;
    RedcIncomplete::call(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    HPBC_CLOCKWORK_ASSERT2(minuend < n && subtrahend < n);

#if defined(HURCHALLA_AVOID_CSELECT)
    T result = static_cast<T>(minuend - subtrahend);
      // T final_result = (minuend < subtrahend) ? result + n : result;
    T mask = static_cast<T>(-static_cast<T>(minuend < subtrahend));
    T final_result = static_cast<T>(result + (mask & n));
#else
# if 0
    T result = static_cast<T>(minuend - subtrahend);
    bool ovf = (minuend < subtrahend);
      // T final_result = ovf ? result + n : result;
    T final_result = hc::conditional_select(ovf, static_cast<T>(result+n), result);
# else
    T final_result = hc::modular_subtraction_prereduced_inputs<T, PTAG>(minuend, subtrahend, n);
# endif
#endif

    HPBC_CLOCKWORK_POSTCONDITION2(final_result < n);
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
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, PTAG());
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }
};


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

# if (HURCHALLA_COMPILER_HAS_UINT128_T())
// specialization for __uint128_t (for x86_64)
template <>
struct RedcStandard<__uint128_t>
{
  using T = __uint128_t;

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowlatencyTag)
  {
    using P = typename safely_promote_unsigned<T>::type;
    // see uint64_t version's comments for explanations
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));
    T mn_hi = ::hurchalla::unsigned_multiply_to_hi_product(m, n);
    HPBC_CLOCKWORK_ASSERT2(mn_hi < n);
    T reg = u_hi + n;

    using std::uint64_t;
    uint64_t reglo = static_cast<uint64_t>(reg);
    uint64_t reghi = static_cast<uint64_t>(reg >> 64);
    uint64_t uhilo = static_cast<uint64_t>(u_hi);
    uint64_t uhihi = static_cast<uint64_t>(u_hi >> 64);
    uint64_t mnhilo = static_cast<uint64_t>(mn_hi);
    uint64_t mnhihi = static_cast<uint64_t>(mn_hi >> 64);
    __asm__ (
        "subq %[mnhilo], %[reglo] \n\t"     /* reg = u_hi + n - mn_hi */
        "sbbq %[mnhihi], %[reghi] \n\t"
        "subq %[mnhilo], %[uhilo] \n\t"     /* t_hi = u_hi - mn_hi */
        "sbbq %[mnhihi], %[uhihi] \n\t"
        "cmovaeq %[uhilo], %[reglo] \n\t"   /* reg = (u_hi >= mn_hi) ? t_hi : reg */
        "cmovaeq %[uhihi], %[reghi] \n\t"
        : [reglo]"+&r"(reglo), [reghi]"+&r"(reghi), [uhilo]"+&r"(uhilo), [uhihi]"+&r"(uhihi)
        : [mnhilo]"r"(mnhilo), [mnhihi]"r"(mnhihi)
        : "cc");
    T result = (static_cast<__uint128_t>(reghi) << 64) | reglo;
    HPBC_CLOCKWORK_ASSERT2(result == DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowlatencyTag()));
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowuopsTag)
  {
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowuopsTag());
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }
};
# endif


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
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));
    T mn_hi = ::hurchalla::unsigned_multiply_to_hi_product(m, n);
    HPBC_CLOCKWORK_ASSERT2(mn_hi < n);
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
    HPBC_CLOCKWORK_ASSERT2(result == DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowlatencyTag()));
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }

  // This function should have: cycles latency 10, fused uops 6
  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowuopsTag)
  {
    // Calling DefaultRedcStandard::call will give us optimal code (we're
    // relying upon modular_subtract_prereduced_inputs() being optimized for low
    // uops - which it is, at least at the time of writing this)
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowuopsTag());
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }
};


// specialization for uint32_t
template <>
struct RedcStandard<std::uint32_t>
{
  using T = std::uint32_t;

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowlatencyTag)
  {
    using P = typename safely_promote_unsigned<T>::type;
    // see uint64_t version's comments for explanations
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));
    T mn_hi = ::hurchalla::unsigned_multiply_to_hi_product(m, n);
    HPBC_CLOCKWORK_ASSERT2(mn_hi < n);
    T reg = u_hi + n;
    T uhi = u_hi;
    __asm__ (
        "subl %[mnhi], %[reg] \n\t"     /* reg = u_hi + n - mn_hi */
        "subl %[mnhi], %[uhi] \n\t"     /* t_hi = u_hi - mn_hi */
        "cmovael %[uhi], %[reg] \n\t"   /* reg = (u_hi >= mn_hi) ? t_hi : reg */
        : [reg]"+&r"(reg), [uhi]"+&r"(uhi)
        : [mnhi]"r"(mn_hi)
        : "cc");
    T result = reg;
    HPBC_CLOCKWORK_ASSERT2(result == DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowlatencyTag()));
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowuopsTag)
  {
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowuopsTag());
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }
};
#endif   // (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) ||
         //  defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) &&
         // defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)







#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) && \
    defined(HURCHALLA_TARGET_ISA_ARM_64) && !defined(_MSC_VER)

// specialization for __uint128_t, if the compiler supports __uint128_t
# if (HURCHALLA_COMPILER_HAS_UINT128_T())
template <>
struct RedcStandard<__uint128_t>
{
  using T = __uint128_t;

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowlatencyTag)
  {
    using P = typename safely_promote_unsigned<T>::type;
    // see uint64_t version's comments for explanations
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));
    T mn_hi = ::hurchalla::unsigned_multiply_to_hi_product(m, n);
    HPBC_CLOCKWORK_ASSERT2(mn_hi < n);
    T reg = u_hi + n;

    using std::uint64_t;
    uint64_t reglo = static_cast<uint64_t>(reg);
    uint64_t reghi = static_cast<uint64_t>(reg >> 64);
    uint64_t uhilo = static_cast<uint64_t>(u_hi);
    uint64_t uhihi = static_cast<uint64_t>(u_hi >> 64);
    uint64_t mnhilo = static_cast<uint64_t>(mn_hi);
    uint64_t mnhihi = static_cast<uint64_t>(mn_hi >> 64);
    __asm__ (
        "subs %[reglo], %[reglo], %[mnhilo] \n\t"       /* reg = u_hi + n - mn_hi */
        "sbcs %[reghi], %[reghi], %[mnhihi] \n\t"
        "subs %[mnhilo], %[uhilo], %[mnhilo] \n\t"      /* res = u_hi - mn_hi */
        "sbcs %[mnhihi], %[uhihi], %[mnhihi] \n\t"
        "csel %[mnhilo], %[reglo], %[mnhilo], lo \n\t"  /* res = (u_hi < mn_hi) ? reg : res */
        "csel %[mnhihi], %[reghi], %[mnhihi], lo \n\t"
        : [reglo]"+&r"(reglo), [reghi]"+&r"(reghi), [mnhilo]"+&r"(mnhilo), [mnhihi]"+&r"(mnhihi)
        : [uhilo]"r"(uhilo), [uhihi]"r"(uhihi)
        : "cc");
    T result = (static_cast<__uint128_t>(mnhihi) << 64) | mnhilo;

    HPBC_CLOCKWORK_ASSERT2(result == DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowlatencyTag()));
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowuopsTag)
  {
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowuopsTag());
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }
};
# endif


// specialization for uint64_t
template <>
struct RedcStandard<std::uint64_t>
{
  using T = std::uint64_t;

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
    HPBC_CLOCKWORK_PRECONDITION2(u_hi < n);
    HPBC_CLOCKWORK_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(inv_n));
    T mn_hi = ::hurchalla::unsigned_multiply_to_hi_product(m, n);
    HPBC_CLOCKWORK_ASSERT2(mn_hi < n);
    T reg = u_hi + n;

    __asm__ (
        "sub %[reg], %[reg], %[mn_hi] \n\t"         /* reg = u_hi + n - mn_hi */
        "subs %[mn_hi], %[u_hi], %[mn_hi] \n\t"     /* res = u_hi - mn_hi */
        "csel %[mn_hi], %[reg], %[mn_hi], lo \n\t"  /* res = (u_hi < mn_hi) ? reg : res */
        : [reg]"+&r"(reg), [mn_hi]"+&r"(mn_hi)
        : [u_hi]"r"(u_hi)
        : "cc");
    T result = mn_hi;

    HPBC_CLOCKWORK_ASSERT2(result == DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowlatencyTag()));
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }

  static HURCHALLA_FORCE_INLINE
  T call(T u_hi, T u_lo, T n, T inv_n, LowuopsTag)
  {
    // Calling DefaultRedcStandard::call will give us optimal code (we're
    // relying upon modular_subtract_prereduced_inputs() being optimized for low
    // uops - which it is, at least at the time of writing this)
    T result = DefaultRedcStandard<T>::call(u_hi, u_lo, n, inv_n, LowuopsTag());
    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }
};

#endif   // (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) ||
         //  defined(HURCHALLA_ALLOW_INLINE_ASM_REDC)) &&
         // defined(HURCHALLA_TARGET_ISA_ARM_64) && !defined(_MSC_VER)





#else


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
    namespace hc = ::hurchalla;

    T minuend, subtrahend;
    RedcIncomplete::call(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    HPBC_CLOCKWORK_ASSERT2(minuend < n && subtrahend < n);
#if 0
    // By RedcIncomplete::call()'s Postcondition #1, we would have:
    T difference = static_cast<T>(minuend - subtrahend);
    bool ovf = (minuend < subtrahend);
        // T result = (ovf) ? static_cast<T>(difference + n) : difference;
    T result = hc::conditional_select((ovf), static_cast<T>(difference + n), difference);
#else
    // The #if section above is just a modular subtraction...
    // The most efficient way to compute it is to call our dedicated function:
    T result = hc::modular_subtraction_prereduced_inputs<T, PTAG>(minuend, subtrahend, n);
#endif

    HPBC_CLOCKWORK_POSTCONDITION2(result < n);
    return result;
  }
};


#endif



}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
