// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTSUB_SQRT_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTSUB_SQRT_RANGE_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// Note: this file is extremely closely related to
// hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_subtraction.h
// However, the allowable input/output ranges differ, which slightly changes the
// arithmetic and necessitates this file.

// The name "sqrt_range" signifies that this function is intended to be used
// with MontySqrtRange.h.

// montsub_sqrt_range::call() requires an unusual input range:
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


// minor note: uses a static member function to disallow ADL.
struct default_montsub_sqrt_range {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
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

    // This implementation is designed for low uop count and low register use.
    // An implementation is possible with expected lower latency and higher uop
    // count and higher register use, but it's not preferred (see older git
    // commits of this file for this alternative).
    T diff = static_cast<T>(a - b);
# if 0
    T result = (a <= b) ? static_cast<T>(diff + n) : diff;
# else
    T result = diff;
    HURCHALLA_CMOV(a <= b, result, static_cast<T>(diff + n));
# endif

    HPBC_POSTCONDITION2(0 < result && result <= n);
    return result;
  }
};


// primary template
template <typename T>
struct montsub_sqrt_range {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T n)
  {
    return default_montsub_sqrt_range::call(a, b, n);
  }
};

#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MONTSUB_SQRT_RANGE)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
template <>
struct montsub_sqrt_range<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t a, std::uint64_t b, std::uint64_t n)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(n > 0);
    // See the discussion in the primary template above regarding the next
    // preconditions, which allow a==0 and/or b==0 so long as we don't have the
    // combination a==0 with b==n.
    HPBC_PRECONDITION2(b <= n);  // 0 <= b is guaranteed by uint64_t
    HPBC_PRECONDITION2(a <= n);  // 0 <= a is guaranteed by uint64_t
    HPBC_PRECONDITION2(!(a == 0 && b == n));

    // Note: the issues and solutions with LEA and RBP/EBP/R13 are the same here
    // as described in impl_modular_subtraction.h
    uint64_t dummy;
    uint64_t result = a;  // we prefer not to overwrite an input (a)
    __asm__ ("subq %[b], %[res] \n\t"            /* res = a - b */
             "leaq (%[res], %[n]), %[tmp] \n\t"  /* tmp = res + n */
             "cmovbeq %[tmp], %[res] \n\t"       /* res = (a<=b) ? tmp : res */

# if defined(__INTEL_COMPILER)
                 : [res]"+&abcdSD"(result), [tmp]"=r"(dummy)
# elif defined(__clang__)    /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
                             /* clang seems to use the first register listed. */
                             /* rcx is probably a good first choice. */
                 : [res]"+&cabdSD"(result), [tmp]"=r"(dummy)
# else
                 : [res]"+&UabcdSD"(result), [tmp]"=r"(dummy)
# endif

# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b]"r"(b), [n]"r"(n)
# else
             : [b]"rm"(b), [n]"r"(n)
# endif

             : "cc");

    HPBC_POSTCONDITION2(0 < result && result <= n);
    HPBC_POSTCONDITION2(result == default_montsub_sqrt_range::call(a, b, n));
    return result;
  }
};
#endif


}} // end namespace

#endif
