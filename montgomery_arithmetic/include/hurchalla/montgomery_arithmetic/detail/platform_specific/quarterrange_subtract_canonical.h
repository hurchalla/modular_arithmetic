// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_QUARTERRANGE_SUBTRACT_CANONICAL_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_QUARTERRANGE_SUBTRACT_CANONICAL_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// ------ quarterrange_subtract_canonical::call(x, y, n) ------
// Intended for use solely by MontyQuarterRange.h
// Returns x-y (mod n).  Note that the return value isn't necessarily canonical.
// y must be canonical (meaning: 0 <= y < n).


// minor note: uses a static member function to disallow ADL.
struct default_quarterrange_subtract_canonical {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T x, T cy, T n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    // We assume this is only called by MontyQuarterRange, which requires
    // that the modulus n < R/4
    constexpr T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                   (ut_numeric_limits<T>::digits - 2));
    HPBC_PRECONDITION2(n < Rdiv4);
    HPBC_PRECONDITION2(x < static_cast<T>(2*n));
    HPBC_PRECONDITION2(cy < n);  // the subtrahend must be canonical

    using S = typename extensible_make_signed<T>::type;
    static_assert(static_cast<S>(static_cast<T>(static_cast<S>(-1))) ==
                  static_cast<S>(-1), "Casting a signed S value to unsigned and"
                               " back again must result in the original value");
    static_assert(static_cast<S>(-1) == ~(static_cast<S>(0)),
                                  "S must use two's complement representation");

    S tmp = static_cast<S>(static_cast<S>(x) - static_cast<S>(cy));
# if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<S>(-1) >> 1) == static_cast<S>(-1),
                      "Arithmetic right shift is required but unavailable");
    // if tmp is negative, create a bit mask of all 1s.  Otherwise all 0s.
    T mask = static_cast<T>(tmp >> ut_numeric_limits<S>::digits);
    T n_masked = static_cast<T>(n & mask);
    T result = static_cast<T>(static_cast<T>(tmp) + n_masked);
# else
    T result = static_cast<T>(static_cast<T>(tmp) + n);
      // result = (tmp>=0) ? tmp : result
    HURCHALLA_CSELECT(result, tmp >= 0, static_cast<T>(tmp), result);
# endif

    HPBC_POSTCONDITION2(result < static_cast<T>(2*n));
    return result;
  }
};


// primary template
template <typename T>
struct quarterrange_subtract_canonical {
  HURCHALLA_FORCE_INLINE static T call(T x, T y, T n)
  {
    return default_quarterrange_subtract_canonical::call(x, y, n);
  }
};


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_QUARTERRANGE_SUBTRACT_CANONICAL)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
template <>
struct quarterrange_subtract_canonical<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t x, std::uint64_t y, std::uint64_t n)
  {
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical
    HPBC_PRECONDITION2(x < static_cast<std::uint64_t>(2*n));

    // We use the "UabcdSD" constraint below so that the LEA instruction doesn't
    // use RBP/EBP or R13 for the base register (which results in slow lea).
    // See impl_modular_subtraction.h for more info.
    std::uint64_t tmp = x;
    std::uint64_t result;
    __asm__ ("subq %[y], %[tmp] \n\t"            /* tmp = x - y */
             "leaq (%[tmp], %[n]), %[res] \n\t"  /* res = tmp + n */
             "cmovaeq %[tmp], %[res] \n\t"       /* res = (x>=y) ? tmp : res */

# if defined(__INTEL_COMPILER)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
# elif defined(__clang__)    /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
                             /* clang seems to use the first register listed. */
                             /* rcx is probably a good first choice. */
             : [tmp]"+&cabdSD"(tmp), [res]"=r"(result)
# else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
# endif

# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [n]"r"(n), [y]"r"(y)
# else
             : [n]"r"(n), [y]"rm"(y)
# endif

             : "cc");

    HPBC_POSTCONDITION2(result < static_cast<std::uint64_t>(2*n));
    HPBC_POSTCONDITION2(result ==
                        default_quarterrange_subtract_canonical::call(x, y, n));
    return result;
  }
};

template <>
struct quarterrange_subtract_canonical<std::uint32_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint32_t call(std::uint32_t x, std::uint32_t y, std::uint32_t n)
  {
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical
    HPBC_PRECONDITION2(x < static_cast<std::uint32_t>(2*n));

    // Regarding the "UabcdSD" constraint, see the uint64_t specialization above
    std::uint32_t tmp = x;
    std::uint32_t result;
    __asm__ ("subl %[y], %[tmp] \n\t"             /* tmp = x - y */
             "leal (%q[tmp], %q[n]), %[res] \n\t" /* res = tmp + n */
             "cmovael %[tmp], %[res] \n\t"        /* res = (x>=y) ? tmp : res */

# if defined(__INTEL_COMPILER)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
# elif defined(__clang__)    /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
                             /* clang seems to use the first register listed. */
                             /* rcx is probably a good first choice. */
             : [tmp]"+&cabdSD"(tmp), [res]"=r"(result)
# else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
# endif

# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [n]"r"(n), [y]"r"(y)
# else
             : [n]"r"(n), [y]"rm"(y)
# endif

             : "cc");

    HPBC_POSTCONDITION2(result < static_cast<std::uint32_t>(2*n));
    HPBC_POSTCONDITION2(result ==
                        default_quarterrange_subtract_canonical::call(x, y, n));
    return result;
  }
};
#endif


}} // end namespace

#endif
