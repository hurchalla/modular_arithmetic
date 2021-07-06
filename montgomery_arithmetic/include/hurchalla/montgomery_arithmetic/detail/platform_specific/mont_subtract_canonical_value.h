// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_SUBTRACT_CANONICAL_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_SUBTRACT_CANONICAL_VALUE_H_INCLUDED


#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <algorithm>

namespace hurchalla { namespace detail {


// mont_subtract_canonical_value()  returns x-y (mod n).
// y must be canonical (meaning: 0 <= y < n).
// The return value is not necessarily canonical, but it is less than or equal
// to max(x, n-1).
template <typename T>
HURCHALLA_FORCE_INLINE T mont_subtract_canonical_value(T x, T y, T n)
{
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical

    T diff = static_cast<T>(x - y);
    T result = static_cast<T>(diff + n);
    // encourage compiler to use conditional move via ternary operator
    result = (x>=y) ? diff : result;

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    return result;
}


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MONT_SUBTRACT_CANONICAL)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// -----------------------------------------------------------------------------
// These non-template functions have first priority for argument matching.
// They have the same descriptions as the generic template version above.
// -----------------------------------------------------------------------------
HURCHALLA_FORCE_INLINE std::uint64_t mont_subtract_canonical_value(
                              std::uint64_t x, std::uint64_t y, std::uint64_t n)
{
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical

    // We use the "UabcdSD" constraint below so that the LEA instruction doesn't
    // use RBP/EBP or R13 for the base register (which results in slow lea).
    // See impl_modular_subtraction.h for more info.
    std::uint64_t tmp = x;
    std::uint64_t result;
    __asm__ ("subq %[y], %[tmp] \n\t"            /* tmp = x - y */
             "leaq (%[tmp], %[n]), %[res] \n\t"  /* res = tmp + n */
             "cmovaeq %[tmp], %[res] \n\t"       /* res = (x>=y) ? tmp : res */
#  if defined(__INTEL_COMPILER) || defined(__clang__)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
#  else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
#  endif
             : [n]"r"(n), [y]"rm"(y)
             : "cc");
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint64_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                         mont_subtract_canonical_value<std::uint64_t>(x, y, n));
    return result;
}

HURCHALLA_FORCE_INLINE std::uint32_t mont_subtract_canonical_value(
                              std::uint32_t x, std::uint32_t y, std::uint32_t n)
{
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical

    // Regarding the "UabcdSD" constraint, see the the uint64_t overload above
    std::uint32_t tmp = x;
    std::uint32_t result;
    __asm__ ("subl %[y], %[tmp] \n\t"             /* tmp = x - y */
             "leal (%q[tmp], %q[n]), %[res] \n\t" /* res = tmp + n */
             "cmovael %[tmp], %[res] \n\t"        /* res = (x>=y) ? tmp : res */
#  if defined(__INTEL_COMPILER) || defined(__clang__)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
#  else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
#  endif
             : [n]"r"(n), [y]"rm"(y)
             : "cc");
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint32_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                         mont_subtract_canonical_value<std::uint32_t>(x, y, n));
    return result;
}
#endif


}} // end namespace

#endif
