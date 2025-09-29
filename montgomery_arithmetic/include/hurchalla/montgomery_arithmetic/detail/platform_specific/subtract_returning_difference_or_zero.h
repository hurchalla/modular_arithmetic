// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_SUBTRACT_RETURNING_DIFFERENCE_OR_ZERO_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_SUBTRACT_RETURNING_DIFFERENCE_OR_ZERO_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// This function in this file is writted specifically to help implement
// MontyFullRange.h's squareSV() function without branches (assuming macros
// are defined to allow assembly here).
// It's unlikely it would be useful for other purposes, but it would be okay
// to do so if it happens to be useful.


// Calculates a - b, placing the subtraction result in 'difference'.
// If a < b, returns difference; otherwise returns 0.
//
// Note: the non-template versions below are always used in C++ instead of this
// template version when there are exact argument matches.  This is good - we
// want to use the asm versions when asm is available.
template <typename T> HURCHALLA_FORCE_INLINE
T subtract_returning_difference_or_zero(T& difference, T a, T b)
{
    namespace hc = ::hurchalla;
    difference = static_cast<T>(a - b);
      // T ret = (a < b) ? difference : static_cast<T>(0);
    T ret = hc::conditional_select((a < b), difference, static_cast<T>(0));
    return ret;
}



// X86_64
// MSVC doesn't support inline asm so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_SUBTRACT_RDZ)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

# if (HURCHALLA_COMPILER_HAS_UINT128_T())
HURCHALLA_FORCE_INLINE __uint128_t
subtract_returning_difference_or_zero(__uint128_t& difference, __uint128_t a, __uint128_t b)
{
    uint64_t diff_lo = static_cast<uint64_t>(a);
    uint64_t diff_hi = static_cast<uint64_t>(a >> 64);
    uint64_t b_lo = static_cast<uint64_t>(b);
    uint64_t b_hi = static_cast<uint64_t>(b >> 64);
    uint64_t ret_lo = 0;
    uint64_t ret_hi = 0;
    __asm__ ("subq %[b_lo], %[diff_lo] \n\t"       /* diff = a - b */
             "sbbq %[b_hi], %[diff_hi] \n\t"
             "cmovbq %[diff_lo], %[ret_lo] \n\t"   /* ret = (a < b) ? diff : 0 */
             "cmovbq %[diff_hi], %[ret_hi] \n\t"
             : [diff_lo]"+&r"(diff_lo), [diff_hi]"+r"(diff_hi),
               [ret_lo]"+r"(ret_lo), [ret_hi]"+r"(ret_hi)
#  if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b_lo]"r"(b_lo), [b_hi]"r"(b_hi)
#  else
             : [b_lo]"rm"(b_lo), [b_hi]"rm"(b_hi)
#  endif
             : "cc");
    difference = (static_cast<__uint128_t>(diff_hi) << 64) | diff_lo;
    __uint128_t ret = (static_cast<__uint128_t>(ret_hi) << 64) | ret_lo;

    HPBC_CLOCKWORK_POSTCONDITION2(difference == a - b);
    HPBC_CLOCKWORK_POSTCONDITION2(ret == (a < b) ? difference : 0);
    return ret;
}
# endif

HURCHALLA_FORCE_INLINE uint64_t
subtract_returning_difference_or_zero(uint64_t& difference, uint64_t a, uint64_t b)
{
    uint64_t diff = a;
    uint64_t ret = 0;
    __asm__ ("subq %[b], %[diff] \n\t"       /* diff = a - b */
             "cmovbq %[diff], %[ret] \n\t"   /* ret = (a < b) ? diff : 0 */
             : [diff]"+r"(diff), [ret]"+r"(ret)
#  if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b]"r"(b)
#  else
             : [b]"rm"(b)
#  endif
             : "cc");
    difference = diff;

    HPBC_CLOCKWORK_POSTCONDITION2(difference == a - b);
    HPBC_CLOCKWORK_POSTCONDITION2(ret == (a < b) ? difference : 0);
    return ret;
}

HURCHALLA_FORCE_INLINE uint32_t
subtract_returning_difference_or_zero(uint32_t& difference, uint32_t a, uint32_t b)
{
    uint32_t diff = a;
    uint32_t ret = 0;
    __asm__ ("subl %[b], %[diff] \n\t"       /* diff = a - b */
             "cmovbl %[diff], %[ret] \n\t"   /* ret = (a < b) ? diff : 0 */
             : [diff]"+r"(diff), [ret]"+r"(ret)
#  if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b]"r"(b)
#  else
             : [b]"rm"(b)
#  endif
             : "cc");
    difference = diff;

    HPBC_CLOCKWORK_POSTCONDITION2(difference == a - b);
    HPBC_CLOCKWORK_POSTCONDITION2(ret == (a < b) ? difference : 0);
    return ret;
}

#endif




// ARM64
// MSVC doesn't support inline asm so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_SUBTRACT_RDZ)) && \
    defined(HURCHALLA_TARGET_ISA_ARM_64) && !defined(_MSC_VER)

# if (HURCHALLA_COMPILER_HAS_UINT128_T())
HURCHALLA_FORCE_INLINE __uint128_t
subtract_returning_difference_or_zero(__uint128_t& difference, __uint128_t a, __uint128_t b)
{
    uint64_t a_lo = static_cast<uint64_t>(a);
    uint64_t a_hi = static_cast<uint64_t>(a >> 64);
    uint64_t b_lo = static_cast<uint64_t>(b);
    uint64_t b_hi = static_cast<uint64_t>(b >> 64);
    uint64_t diff_lo, diff_hi, ret_lo, ret_hi;
    __asm__ ("subs %[diff_lo], %[a_lo], %[b_lo] \n\t"      /* diff = a - b */
             "sbcs %[diff_hi], %[a_hi], %[b_hi] \n\t"
             "csel %[ret_lo], %[diff_lo], xzr, lo \n\t"    /* ret = (a < b) ? diff : 0 */
             "csel %[ret_hi], %[diff_hi], xzr, lo \n\t"
             : [diff_lo]"=&r"(diff_lo), [diff_hi]"=r"(diff_hi),
               [ret_lo]"=r"(ret_lo), [ret_hi]"=r"(ret_hi)
             : [a_lo]"r"(a_lo), [a_hi]"r"(a_hi), [b_lo]"r"(b_lo), [b_hi]"r"(b_hi)
             : "cc");
    difference = (static_cast<__uint128_t>(diff_hi) << 64) | diff_lo;
    __uint128_t ret = (static_cast<__uint128_t>(ret_hi) << 64) | ret_lo;

    HPBC_CLOCKWORK_POSTCONDITION2(difference == a - b);
    HPBC_CLOCKWORK_POSTCONDITION2(ret == (a < b) ? difference : 0);
    return ret;
}
# endif

HURCHALLA_FORCE_INLINE uint64_t
subtract_returning_difference_or_zero(uint64_t& difference, uint64_t a, uint64_t b)
{
    uint64_t diff, ret;
    __asm__ ("subs %[diff], %[a], %[b] \n\t"        /* diff = a - b */
             "csel %[ret], %[diff], xzr, lo \n\t"   /* ret = (a < b) ? diff : 0 */
             : [diff]"=r"(diff), [ret]"=r"(ret)
             : [a]"r"(a), [b]"r"(b)
             : "cc");
    difference = diff;

    HPBC_CLOCKWORK_POSTCONDITION2(difference == a - b);
    HPBC_CLOCKWORK_POSTCONDITION2(ret == (a < b) ? difference : 0);
    return ret;
}

HURCHALLA_FORCE_INLINE uint32_t
subtract_returning_difference_or_zero(uint32_t& difference, uint32_t a, uint32_t b)
{
    uint64_t a64 = a;
    uint64_t b64 = b;
    uint64_t diff64;
    uint64_t ret64 = subtract_returning_difference_or_zero(diff64, a64, b64);
    difference = static_cast<uint32_t>(diff64);
    uint32_t ret = static_cast<uint32_t>(ret64);

    HPBC_CLOCKWORK_POSTCONDITION2(difference == a - b);
    HPBC_CLOCKWORK_POSTCONDITION2(ret == (a < b) ? difference : 0);
    return ret;
}

#endif


}} // end namespace

#endif
