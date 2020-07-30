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

// The name "sqrt_range" signifies that this function is intended to be used
// with MontySqrtRange.h.

// Our function montsub_sqrt_range() requires an unusual input range:
// 0 <= a <= n,  and  0 <= b <= n  (so long as !(a==0 && b==n)).
// The output return value range will be  0 < returnValue <= n.
// Obviously neither the inputs nor outputs necessarily belong to the minimal
// residue class modulo n, since they are allowed to equal n.
// These preconditions and postconditions originate from MontySqrtRange.h,
// although the preconditions here are relaxed slightly from MontySqrtRange.
// They allow this function to be used seemlessly by MontySqrtRange, since
// MontySqrtRange will always provide inputs that respect our preconditions,
// and our postconditions ensure we will always provide valid values for
// MontySqrtRange.


template <typename T>
HURCHALLA_FORCE_INLINE T montsub_sqrt_range(T a, T b, T n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(n > 0);
    // MontySqrtRange uses input/output values that satisfy 0 < value <= n, but
    // we can relax the precondition here to allow zero values for a or b, even
    // though MontySqrtRange won't send those values.  This will be useful in
    // some circumstances, though we need to make sure we'll still be able to
    // correctly satisfy our postcondition of  0 < result <= n.
    // To accomplish this we just need to make sure our newly relaxed
    // preconditions disallow a==0 with b==n, since that is the one and only
    // combination that is problematic.
    HPBC_PRECONDITION2(0 <= b && b <= n);
    HPBC_PRECONDITION2(0 <= a && a <= n);
    HPBC_PRECONDITION2(!(a == 0 && b == n));

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
HURCHALLA_FORCE_INLINE std::uint64_t montsub_sqrt_range(std::uint64_t a,
                                               std::uint64_t b, std::uint64_t n)
{
    using std::uint64_t;
    HPBC_PRECONDITION2(n > 0);
    // See the discussion in the above template function regarding the next
    // preconditions, which allow a==0 and/or b==0 so long as we don't have the
    // combination a==0 with b==n.
    HPBC_PRECONDITION2(b <= n);  // 0 <= b is guaranteed by uint64_t
    HPBC_PRECONDITION2(a <= n);  // 0 <= a is guaranteed by uint64_t
    HPBC_PRECONDITION2(!(a == 0 && b == n));

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
    HPBC_POSTCONDITION2(result == montsub_sqrt_range<uint64_t>(a, b, n));
    return result;
}
#endif


}} // end namespace

#endif
