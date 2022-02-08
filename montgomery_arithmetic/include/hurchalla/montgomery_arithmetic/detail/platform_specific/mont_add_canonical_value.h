// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_ADD_CANONICAL_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_ADD_CANONICAL_VALUE_H_INCLUDED


#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <algorithm>

namespace hurchalla { namespace detail {


// mont_add_canonical_value::call()  returns x+y (mod n).
// y must be canonical (meaning: 0 <= y < n).
// The return value is not necessarily canonical, but it is less than or equal
// to max(x, n-1).


// minor note: uses a static member function to disallow ADL.
struct default_mont_add_canonical_value {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T x, T y, T n)
  {
    HPBC_PRECONDITION2(y < n);  // the second addend must be canonical

    // Naively, we would like to set  result = (x+y >= n) ? (x+y-n) : x+y.
    // But x+y possibly could overflow, so instead we will use the equivalent
    // conditional (x >= n-y).  This is safe since due to the precondition y<n,
    // n-y will never overflow.  So we set
    // result = (x >= n-y) ? (x-(n-y)) : x+y
    T tmp = static_cast<T>(n - y);
    T sum = static_cast<T>(x + y);
    T tmp2 = static_cast<T>(x - tmp);
#if 0
    // encourage compiler to use conditional move via ternary operator
    T result = (x>=tmp) ? tmp2 : sum;
#else
    T result = sum;
    HURCHALLA_CMOV(x >= tmp, result, tmp2);
#endif

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    return result;
  }
};


// primary template
template <typename T>
struct mont_add_canonical_value {
  HURCHALLA_FORCE_INLINE static T call(T x, T y, T n)
  {
    return default_mont_add_canonical_value::call(x, y, n);
  }
};


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MONT_ADD_CANONICAL)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
template <>
struct mont_add_canonical_value<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t x, std::uint64_t y, std::uint64_t n)
  {
    HPBC_PRECONDITION2(y < n);  // the second addend must be canonical
    std::uint64_t tmp = static_cast<std::uint64_t>(n - y);
    std::uint64_t sum = static_cast<std::uint64_t>(x + y);
    std::uint64_t tmp2 = x;
    __asm__ ("subq %[tmp], %[tmp2] \n\t"     /* tmp2 = x - tmp */
             "cmovaeq %[tmp2], %[sum] \n\t"  /* sum = (x>=tmp) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
# if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmp]"r"(tmp)
# else
             : [tmp]"rm"(tmp)
# endif
             : "cc");
    std::uint64_t result = sum;
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint64_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                               default_mont_add_canonical_value::call(x, y, n));
    return result;
  }
};

template <>
struct mont_add_canonical_value<std::uint32_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint32_t call(std::uint32_t x, std::uint32_t y, std::uint32_t n)
  {
    HPBC_PRECONDITION2(y < n);
    std::uint32_t tmp = static_cast<std::uint32_t>(n - y);
    std::uint32_t sum = static_cast<std::uint32_t>(x + y);
    std::uint32_t tmp2 = x;
    __asm__ ("subl %[tmp], %[tmp2] \n\t"      /* tmp2 = x - tmp */
             "cmovael %[tmp2], %[sum] \n\t"   /* sum = (x>=tmp) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
# if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmp]"r"(tmp)
# else
             : [tmp]"rm"(tmp)
# endif
             : "cc");
    std::uint32_t result = sum;
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint32_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                               default_mont_add_canonical_value::call(x, y, n));
    return result;
  }
};
#endif


}} // end namespace

#endif
