// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_ADD_CANONICAL_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_ADD_CANONICAL_VALUE_H_INCLUDED


#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <algorithm>

namespace hurchalla { namespace detail {


// mont_add_canonical_value()  returns x+y (mod n).
// y must be canonical (meaning: 0 <= y < n).
// The return value is not necessarily canonical, but it is less than or equal
// to max(x, n-1).
template <typename T>
HURCHALLA_FORCE_INLINE T mont_add_canonical_value(T x, T y, T n)
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
    // encourage compiler to use conditional move via ternary operator
    T result = (x>=tmp) ? tmp2 : sum;

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    return result;
}


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MONT_ADD_CANONICAL)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// -----------------------------------------------------------------------------
// These non-template functions have first priority for argument matching.
// They have the same descriptions as the generic template version above.
// -----------------------------------------------------------------------------
HURCHALLA_FORCE_INLINE std::uint64_t mont_add_canonical_value(std::uint64_t x,
                                               std::uint64_t y, std::uint64_t n)
{
    HPBC_PRECONDITION2(y < n);  // the second addend must be canonical
    std::uint64_t tmp = static_cast<std::uint64_t>(n - y);
    std::uint64_t sum = static_cast<std::uint64_t>(x + y);
    std::uint64_t tmp2 = x;
    __asm__ ("subq %[tmp], %[tmp2] \n\t"     /* tmp2 = x - tmp */
             "cmovaeq %[tmp2], %[sum] \n\t"  /* sum = (x>=tmp) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
             : [tmp]"rm"(tmp)
             : "cc");
    std::uint64_t result = sum;
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint64_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                              mont_add_canonical_value<std::uint64_t>(x, y, n));
    return result;
}

HURCHALLA_FORCE_INLINE std::uint32_t mont_add_canonical_value(std::uint32_t x,
                                               std::uint32_t y, std::uint32_t n)
{
    HPBC_PRECONDITION2(y < n);
    std::uint32_t tmp = static_cast<std::uint32_t>(n - y);
    std::uint32_t sum = static_cast<std::uint32_t>(x + y);
    std::uint32_t tmp2 = x;
    __asm__ ("subl %[tmp], %[tmp2] \n\t"      /* tmp2 = x - tmp */
             "cmovael %[tmp2], %[sum] \n\t"   /* sum = (x>=tmp) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
             : [tmp]"rm"(tmp)
             : "cc");
    std::uint32_t result = sum;
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint32_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                              mont_add_canonical_value<std::uint32_t>(x, y, n));
    return result;
}
#endif


}} // end namespace

#endif
