// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_HALFRANGE_GET_CANONICAL_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_HALFRANGE_GET_CANONICAL_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <algorithm>

namespace hurchalla { namespace detail {


// ------ halfrange_get_canonical::call(x, n) ------
// Intended for use solely by MontyHalfRange.h
// Returns  result = x (mod n),  with 0 <= result < n.


// minor note: uses a static member function to disallow ADL.
struct default_halfrange_get_canonical {
  template <typename S>
  HURCHALLA_FORCE_INLINE static S call(S x, S n)
  {
    // We assume this function is only called by MontyHalfRange
    static_assert(ut_numeric_limits<S>::is_integer, "");
    static_assert(ut_numeric_limits<S>::is_signed, "");
    using T = typename extensible_make_unsigned<S>::type;
    static_assert(static_cast<S>(-1) == ~(static_cast<S>(0)),
                                  "S must use two's complement representation");
    static_assert(static_cast<S>(static_cast<T>(static_cast<S>(-1))) ==
                  static_cast<S>(-1), "Casting a signed S value to unsigned and"
                               " back again must result in the original value");
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(-n <= x && x < n);

    T tx = static_cast<T>(x);
    T tn = static_cast<T>(n);
#if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<S>(-1) >> 1) == static_cast<S>(-1),
                          "Arithmetic right shift is required but unavailable");
    // Use arithmetic right shift of sign bit to create mask of all 1s or 0s
    T mask = static_cast<T>(x >> ut_numeric_limits<S>::digits);
    T n_masked = static_cast<T>(tn & mask);
    T tc = static_cast<T>(tx + n_masked);
#else
    T tc = static_cast<T>(tx + tn);
#  if defined(__clang__)
    // Since we know  -n <= x < n, the sum for tc above will have carried
    // (wrapped-around) if x<0.  Thus testing (tc < tn) is equivalent to a test
    // for (x < 0).  And consequently, testing (tc >= tn) is equivalent to a
    // test for (x >= 0).  It is in theory the better test, because it lets the
    // compiler use the flags already set by (tx+tn), if the compiler is smart.
    // Only clang appears to be smart enough to do this though; for x64 and
    // arm32/64: gcc, msvc, icc produce better asm with the #else.
    tc = conditional_select(tc >= tn, tx, tc);  // tc = (tc>=n) ? tx : tc
#  else
    tc = conditional_select(x >= 0, tx, tc);  // tc = (x>=0) ? tx : tc
#  endif
#endif

    HPBC_POSTCONDITION2(tc < tn);
    return static_cast<S>(tc);
  }
};


// primary template
template <typename S>
struct halfrange_get_canonical {
  HURCHALLA_FORCE_INLINE static S call(S x, S n)
  {
    return default_halfrange_get_canonical::call(x, n);
  }
};


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_HALFRANGE_GET_CANONICAL)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
template <>
struct halfrange_get_canonical<std::int64_t> {
  HURCHALLA_FORCE_INLINE
  static std::int64_t call(std::int64_t x, std::int64_t n)
  {
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(-n <= x && x < n);

    std::int64_t tmp = x;
    __asm__ ("addq %[n], %[tmp] \n\t"       /* tmp = x + n */
             "cmovaeq %[x], %[tmp] \n\t"    /* tmp = (tmp>=n) ? x : tmp */
             : [tmp]"+&r"(tmp)
# if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [n]"r"(n), [x]"r"(x)
# else
             : [n]"rm"(n), [x]"rm"(x)
# endif
             : "cc");
    std::int64_t result = tmp;
    HPBC_POSTCONDITION2(0 <= result && result < n);
    HPBC_POSTCONDITION2(result == default_halfrange_get_canonical::call(x,n));
    return result;
  }
};

template <>
struct halfrange_get_canonical<std::int32_t> {
  HURCHALLA_FORCE_INLINE
  static std::int32_t call(std::int32_t x, std::int32_t n)
  {
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(-n <= x && x < n);

    std::int32_t tmp = x;
    __asm__ ("addl %[n], %[tmp] \n\t"       /* tmp = x + n */
             "cmovael %[x], %[tmp] \n\t"    /* tmp = (tmp>=n) ? x : tmp */
             : [tmp]"+&r"(tmp)
# if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [n]"r"(n), [x]"r"(x)
# else
             : [n]"rm"(n), [x]"rm"(x)
# endif
             : "cc");
    std::int32_t result = tmp;
    HPBC_POSTCONDITION2(0 <= result && result < n);
    HPBC_POSTCONDITION2(result == default_halfrange_get_canonical::call(x,n));
    return result;
  }
};
#endif


}} // end namespace

#endif
