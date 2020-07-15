// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTSUB_SQRT_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTSUB_SQRT_RANGE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


// Note: this file is extremely closely related to
// hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_subtraction.h
// However, the allowable input/output ranges differ, which slightly changes the
// arithmetic and necessitates this file.

// The name "sqrt_range" signifies that modulus 'n' must be less than sqrt(R),
// where R = 2^(ma_numeric_limits<T>::digits).  For example, if T is uint64_t
// then R = 2^64 and sqrt(R) == 2^32, and thus we would require  n < 2^32.


template <typename T>
HURCHALLA_FORCE_INLINE T montsub_sqrt_range(T a, T b, T n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(0 < a && a <= n);
    HPBC_PRECONDITION2(0 < b && b <= n);

    // We want essentially-  result = (a-b <= 0) ? a-b+n : a-b
    //    But (a-b) overflows whenever b>a, so instead of testing if (a-b <= 0),
    //   we test the alternative predicate (a <= b).  This gives us our desired
    //   result without any problem of overflow.  So we can and should use:
    //   result = (a <= b) ? a-b+n : a-b

    // we calculate diff here to encourage the compiler to hoist it out of a
    // loop (assuming this function is inlined inside a loop).
    T diff = static_cast<T>(n - b);

    T tmp = static_cast<T>(a + diff);
    T tmp2 = static_cast<T>(a - b);
    T result = (a <= b) ? tmp : tmp2;

    HPBC_POSTCONDITION2(0 < result && result <= n);
    return result;
}


// -------- PLATFORM SPECIFIC nontemplate overloads ----------

// Note: a nontemplate function overload gets first priority for being called
// (see http://www.gotw.ca/publications/mill17.htm ), when both the nontemplate
// function and the generic template function match the caller's provided
// argument type(s).

#if defined(HURCHALLA_ALLOW_INLINE_ASM_MODSUB) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// This function is an asm version of the template montsub_sqrt_range()
HURCHALLA_FORCE_INLINE uint64_t montsub_sqrt_range(uint64_t a, uint64_t b,
                                                                     uint64_t n)
{
    HPBC_PRECONDITION2(n > 0);
    HPBC_PRECONDITION2(0 < a && a <= n);
    HPBC_PRECONDITION2(0 < b && b <= n);

    // By calculating diff outside of the __asm__, we allow the compiler to loop
    // hoist diff, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint64_t diff = n - b;
    uint64_t tmp = a + diff;

    uint64_t result;
    __asm__ ("subq %[b], %0 \n\t"       /* result = a - b */
             "cmovbeq %[tmp], %0 \n\t"  /* result = (a <= b) ? tmp : result */
             : "=&r"(result)
             : "0"(a), [b]"r"(b), [tmp]"r"(tmp)
             : "cc");

    HPBC_POSTCONDITION2(0 < result && result <= n);
    return result;
}
#endif


}} // end namespace

#endif
