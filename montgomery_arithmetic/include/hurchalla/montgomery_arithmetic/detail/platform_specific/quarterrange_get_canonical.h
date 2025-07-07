// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_QUARTERRANGE_GET_CANONICAL_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_QUARTERRANGE_GET_CANONICAL_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <algorithm>

namespace hurchalla { namespace detail {


// ------ quarterrange_get_canonical::call(x, n) ------
// Intended for use solely by MontyQuarterRange.h
// Returns  result = x (mod n),  with 0 <= result < n.


// Fyi: we use structs with static member functions to disallow ADL and to make
// specializations simple and easy.

struct default_quarterrange_get_canonical {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T x, T n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    // We assume this is only called by MontyQuarterRange, which requires
    // that the modulus n < R/4
    constexpr T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                   (ut_numeric_limits<T>::digits - 2));
    HPBC_PRECONDITION2(0 <= n && n < Rdiv4);
    HPBC_PRECONDITION2(0 <= x && x < static_cast<T>(2*n));

#if 0
   // this should be correct, but the #else is preferred for performance
    T result = static_cast<T>(x - n);
       // set  result = (x<n) ? x : result
    result = ::hurchalla::conditional_select(x<n, x, result);
#else
    using S = typename extensible_make_signed<T>::type;
    static_assert(static_cast<S>(-1) == ~(static_cast<S>(0)),
                                  "S must use two's complement representation");
    static_assert(static_cast<S>(static_cast<T>(static_cast<S>(-1))) ==
                  static_cast<S>(-1), "Casting a signed S value to unsigned and"
                               " back again must result in the original value");
    S tmp = static_cast<S>(static_cast<S>(x) - static_cast<S>(n));
    // tmp has range [-n, n)
# if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<S>(-1) >> 1) == static_cast<S>(-1),
                          "Arithmetic right shift is required but unavailable");
    // if tmp is negative, create a bit mask of all 1s.  Otherwise all 0s.
    T mask = static_cast<T>(tmp >> ut_numeric_limits<S>::digits);
    T n_masked = static_cast<T>(n & mask);
    T result = static_cast<T>(static_cast<T>(tmp) + n_masked);
# else
       // tmp = (tmp<0) ? x : tmp
    tmp = ::hurchalla::conditional_select(tmp<0, static_cast<S>(x), tmp);
    HPBC_ASSERT2(0 <= tmp && tmp < static_cast<S>(n));
    T result = static_cast<T>(tmp);
# endif
#endif
    HPBC_POSTCONDITION2(result < n);
    return result;
  }
};


// primary template
template <typename T>
struct quarterrange_get_canonical {
  HURCHALLA_FORCE_INLINE static T call(T x, T n)
  {
    return default_quarterrange_get_canonical::call(x, n);
  }
};


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_QUARTERRANGE_GET_CANONICAL)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

#ifdef HURCHALLA_ENABLE_INLINE_ASM_128_BIT
template <>
struct quarterrange_get_canonical<__uint128_t> {
  HURCHALLA_FORCE_INLINE
  static __uint128_t call(__uint128_t x, __uint128_t n)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(x < static_cast<__uint128_t>(2*n));

    uint64_t xlo = static_cast<uint64_t>(x);
    uint64_t xhi = static_cast<uint64_t>(x >> 64);
    uint64_t xlo2 = xlo;
    uint64_t xhi2 = xhi;
    uint64_t nlo = static_cast<uint64_t>(n);
    uint64_t nhi = static_cast<uint64_t>(n >> 64);
    __asm__ ("subq %[nlo], %[xlo] \n\t"       /* res = x - n */
             "sbbq %[nhi], %[xhi] \n\t"
             "cmovbq %[xlo2], %[xlo] \n\t"    /* res = (x<n) ? x : res */
             "cmovbq %[xhi2], %[xhi] \n\t"
             : [xlo]"+&r"(xlo), [xhi]"+&r"(xhi)
# if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [nlo]"r"(nlo), [nhi]"r"(nhi), [xlo2]"r"(xlo2), [xhi2]"r"(xhi2)
# else
             : [nlo]"rm"(nlo), [nhi]"rm"(nhi), [xlo2]"rm"(xlo2), [xhi2]"rm"(xhi2)
# endif
             : "cc");
    __uint128_t result = (static_cast<__uint128_t>(xhi) << 64) | xlo;

    HPBC_POSTCONDITION2(result < n);
    HPBC_POSTCONDITION2(result== default_quarterrange_get_canonical::call(x,n));
    return result;
  }
};
#endif

template <>
struct quarterrange_get_canonical<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t x, std::uint64_t n)
  {
    HPBC_PRECONDITION2(x < static_cast<std::uint64_t>(2*n));

    std::uint64_t res = x;
    std::uint64_t tmp = x;
    __asm__ ("subq %[n], %[res] \n\t"       /* res = x - n */
             "cmovbq %[tmp], %[res] \n\t"   /* res = (x<n) ? x : res */
             : [res]"+&r"(res)
# if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [n]"r"(n), [tmp]"r"(tmp)
# else
             : [n]"rm"(n), [tmp]"rm"(tmp)
# endif
             : "cc");
    std::uint64_t result = res;

    HPBC_POSTCONDITION2(result < n);
    HPBC_POSTCONDITION2(result== default_quarterrange_get_canonical::call(x,n));
    return result;
  }
};

template <>
struct quarterrange_get_canonical<std::uint32_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint32_t call(std::uint32_t x, std::uint32_t n)
  {
    HPBC_PRECONDITION2(x < static_cast<std::uint32_t>(2*n));

    std::uint32_t res = x;
    std::uint32_t tmp = x;
    __asm__ ("subl %[n], %[res] \n\t"       /* res = x - n */
             "cmovbl %[tmp], %[res] \n\t"   /* res = (x<n) ? x : res */
             : [res]"+&r"(res)
# if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [n]"r"(n), [tmp]"r"(tmp)
# else
             : [n]"rm"(n), [tmp]"rm"(tmp)
# endif
             : "cc");
    std::uint32_t result = res;

    HPBC_POSTCONDITION2(result < n);
    HPBC_POSTCONDITION2(result== default_quarterrange_get_canonical::call(x,n));
    return result;
  }
};
#endif


}} // end namespace

#endif
