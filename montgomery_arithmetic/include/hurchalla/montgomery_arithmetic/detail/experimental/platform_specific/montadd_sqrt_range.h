// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTADD_SQRT_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTADD_SQRT_RANGE_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// Note: this file is extremely closely related to
// hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_addition.h
// However, the allowable input/output ranges differ, which slightly changes the
// arithmetic and necessitates this file.

// The name "sqrt_range" signifies that this function is intended to be used
// with MontySqrtRange.h.

// montadd_sqrt_range::call() requires/allows an unusual input range:  we allow
// 0 < a <= n,  and  0 < b <= n.
// Similarly, the output return value range will be  0 < returnValue <= n.
// Obviously neither inputs nor outputs necessarily belong to the minimal
// residue class modulo n, since they are allowed to equal n (and not 0).
// These preconditions and postconditions originate from MontySqrtRange.h.
// They allow this function to be used seemlessly by MontySqrtRange, since
// MontySqrtRange will always provide inputs that respect our preconditions,
// and our postconditions ensure we will always provide valid values for
// MontySqrtRange.

// For discussion purposes, let R = 1<<(ut_numeric_limits<T>::digits).  For
// example if T is uint64_t, then R = 1<<64.


// minor note: uses a static member function to disallow ADL.
struct default_montadd_sqrt_range {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(0 < a && a <= n);
    HPBC_PRECONDITION2(0 < b && b <= n);

    // We want essentially-  result = (a+b <= n) ? a+b : a+b-n
    //   We don't need to worry about overflow on (a+b) because a<n and b<n
    //   and n<sqrt(R).  However, for consistency with impl_modular_addition.h,
    //   we will test the alternative predicate (a <= n-b), which will give us
    //   our desired result.  This predicate also has the advantage that (n-b)
    //   might be potentially loop hoisted by the compiler, if this function is
    //   inlined into a loop (n and b might be unmodified by the loop, whereas
    //   'a' probably will probably change on each loop iteration).
    //   So we will use:
    //   result = (a <= n-b) ? a+b : a+b-n
    T tmp = static_cast<T>(n - b);
    T sum = static_cast<T>(a + b);
    T tmp2 = static_cast<T>(a - tmp);
    T result;  // set result = (a <= tmp) ? sum : tmp2
    HURCHALLA_CSELECT(result, a <= tmp, sum, tmp2);

    HPBC_POSTCONDITION2(0 < result && result <= n);
    return result;
  }
};






// primary template
template <typename T>
struct montadd_sqrt_range {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T n)
  {
    return default_montadd_sqrt_range::call(a, b, n);
  }
};

#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MONTADD_SQRT_RANGE)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
template <>
struct montadd_sqrt_range<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t a, std::uint64_t b, std::uint64_t n)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(0 < a && a <= n);
    HPBC_PRECONDITION2(0 < b && b <= n);

    // By calculating tmp outside of the __asm__, we allow the compiler to
    // potentially loop hoist tmp, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint64_t tmp = n - b;
    uint64_t sum = a + b;
    uint64_t tmp2 = a;  // we prefer not to overwrite an input (a)
    __asm__ ("subq %[tmp], %[tmp2] \n\t"     /* tmp2 = a - tmp */
             "cmovbeq %[sum], %[tmp2] \n\t"  /* tmp2 = (a<=tmp) ? sum : tmp2 */
             : [tmp2]"+&r"(tmp2)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmp]"r"(tmp), [sum]"r"(sum)
# else
             : [tmp]"rm"(tmp), [sum]"r"(sum)
# endif
             : "cc");
    uint64_t result = tmp2;

    HPBC_POSTCONDITION2(0 < result && result <= n);
    HPBC_POSTCONDITION2(result == default_montadd_sqrt_range::call(a, b, n));
    return result;
  }
};
#endif


}} // end namespace

#endif
