// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTADD_SQRT_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTADD_SQRT_RANGE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


// Note: this file is extremely closely related to
// hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_addition.h
// However, the allowable input/output ranges differ, which slightly changes the
// arithmetic and necessitates this file.

// The name "sqrt_range" signifies that this function is intended to be used
// with MontySqrtRange.h.

// montadd_sqrt_range() requires/allows an unusual input range:  we allow
// 0 < a <= n,  and  0 < b <= n.
// Similarly, the output return value range will be  0 < returnValue <= n.
// Obviously neither inputs nor outputs necessarily belong to the minimal
// residue class modulo n, since they are allowed to equal n (and not 0).
// These preconditions and postconditions originate from MontySqrtRange.h.
// They allow this function to be used seemlessly by MontySqrtRange, since
// MontySqrtRange will always provide inputs that respect our preconditions,
// and our postconditions ensure we will always provide valid values for
// MontySqrtRange.


template <typename T>
HURCHALLA_FORCE_INLINE T montadd_sqrt_range(T a, T b, T n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(0 < a && a <= n);
    HPBC_PRECONDITION2(0 < b && b <= n);

    // We want essentially-  result = (a+b <= n) ? a+b : a+b-n
    //   However, for consistency with impl_modular_addition.h and potentially
    //   to avoid overflow on (a+b), we will instead test the alternative
    //   predicate (a <= n-b), which will give us our desired result.  This
    //   predicate also has the advantage that (n-b) might be potentially loop
    //   hoisted by the compiler, if this function is inlined into a loop (n and
    //   b might be unmodified by the loop, whereas 'a' probably will probably
    //   change on each loop iteration).
    //   So we will use:
    //   result = (a <= n-b) ? a+b : a+b-n
    T tmp = static_cast<T>(n - b);
    T sum = static_cast<T>(a + b);
    T tmp2 = static_cast<T>(a - tmp);
    T result = (a <= tmp) ? sum : tmp2;

    HPBC_POSTCONDITION2(0 < result && result <= n);
    return result;
}


// -------- PLATFORM SPECIFIC nontemplate overloads ----------

// Note: a nontemplate function overload gets first priority for being called
// (see http://www.gotw.ca/publications/mill17.htm ), when both the nontemplate
// function and the generic template function match the caller's provided
// argument type(s).


#if defined(HURCHALLA_ALLOW_INLINE_ASM_MODADD) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// This function is an asm version of the template montadd_sqrt_range()
HURCHALLA_FORCE_INLINE std::uint64_t montadd_sqrt_range(std::uint64_t a,
                                               std::uint64_t b, std::uint64_t n)
{
    using std::uint64_t;
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(0 < a && a <= n);
    HPBC_PRECONDITION2(0 < b && b <= n);

    // By calculating tmp outside of the __asm__, we allow the compiler to loop
    // hoist tmp, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint64_t tmp = n - b;
    uint64_t sum = a + b;
    uint64_t result;
    __asm__ ("subq %[tmp], %0 \n\t"     /* tmp2 = a - tmp */
             "cmovbeq %[sum], %0 \n\t"  /* result = (a<=tmp) ? sum : tmp2 */
             : "=&r"(result)
             : "0"(a), [tmp]"r"(tmp), [sum]"r"(sum)
             : "cc");

    HPBC_POSTCONDITION2(0 < result && result <= n);
    HPBC_POSTCONDITION2(result == montadd_sqrt_range<uint64_t>(a, b, n));
    return result;
}
#endif


}} // end namespace

#endif
