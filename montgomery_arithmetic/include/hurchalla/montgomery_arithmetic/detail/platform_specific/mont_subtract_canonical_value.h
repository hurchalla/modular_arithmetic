// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_SUBTRACT_CANONICAL_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_SUBTRACT_CANONICAL_VALUE_H_INCLUDED


#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <algorithm>

namespace hurchalla { namespace detail {


// mont_subtract_canonical_value::call()  returns x-y (mod n).
// y must be canonical (meaning: 0 <= y < n).
// The return value is not necessarily canonical, but it is less than or equal
// to max(x, n-1).


// minor note: uses a static member function to disallow ADL.
struct default_mont_subtract_canonical_value {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T x, T y, T n)
  {
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical

    T diff = static_cast<T>(x - y);
    T result = static_cast<T>(diff + n);
# if 0
    // encourage compiler to use conditional move via ternary operator
    result = (x>=y) ? diff : result;
# else
    HURCHALLA_CMOV(x>=y, result, diff);
# endif

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    return result;
  }
};


// primary template
template <typename T>
struct mont_subtract_canonical_value {
  HURCHALLA_FORCE_INLINE static T call(T x, T y, T n)
  {
    return default_mont_subtract_canonical_value::call(x, y, n);
  }
};


#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MONT_SUBTRACT_CANONICAL)) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
template <>
struct mont_subtract_canonical_value<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t x, std::uint64_t y, std::uint64_t n)
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

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint64_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                          default_mont_subtract_canonical_value::call(x, y, n));
    return result;
  }
};

template <>
struct mont_subtract_canonical_value<std::uint32_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint32_t call(std::uint32_t x, std::uint32_t y, std::uint32_t n)
  {
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical

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

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<std::uint32_t>(n-1)));
    HPBC_POSTCONDITION2(result ==
                          default_mont_subtract_canonical_value::call(x, y, n));
    return result;
  }
};
#endif


}} // end namespace

#endif
