// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_HALF_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_HALF_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace montgomery_arithmetic {


// The name "half_range" signifies that modulus 'n' must be less than R/2, where
// R = 2^(ma_numeric_limits<T>::digits).  For example, if T is uint64_t then
// R = 2^64 and R/2 == 2^63, and thus we would require  n < 2^63.


template <typename T>
HURCHALLA_FORCE_INLINE T montmul_half_range(T x, T y, T n, T neg_inv_n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(n < Rdiv2);
    }
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);

    // x<n with y<n  will satisfy montmul_non_minimized's precondition
    // requirement that x*y < n*R.
    bool ovf;
    T prod = montmul_non_minimized(ovf, x, y, n, neg_inv_n);

    // Since we know that the modulus  n < R/2, the
    // montmul_non_minimized()'s postconditions guarantee ovf == false.
    HPBC_ASSERT2(ovf == false);
    // Since ovf == false, montmul_non_minimized()'s postconditions guarantee
    T minimized_result = (prod >= n) ? static_cast<T>(prod-n) : prod;

    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
}


// -------- PLATFORM SPECIFIC nontemplate overloads ----------

// Note: a nontemplate function overload gets first priority for being called
// (see http://www.gotw.ca/publications/mill17.htm ), when both the nontemplate
// function and the generic template function match the caller's provided
// argument type(s).


#if defined(HURCHALLA_ALLOW_INLINE_ASM_MONTMUL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// This function is an asm version of the template function montmul_half_range()
HURCHALLA_FORCE_INLINE std::uint64_t montmul_half_range(std::uint64_t x,
                      std::uint64_t y, std::uint64_t n, std::uint64_t neg_inv_n)
{
    using std::uint64_t;
    HPBC_PRECONDITION2(n < (static_cast<uint64_t>(1) << 63));
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);

    // Performing this first multiply outside the asm block lets the compiler
    // potentially optimize it.  E.g. if it can see x == y, or x == 1 or y == 1.
    __uint128_t u = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(y);
    uint64_t u_hi = static_cast<uint64_t>(u >> 64);   // rdx
    uint64_t u_lo = static_cast<uint64_t>(u);         // rax

    // montgomery REDC portion.
    // See montmul_half_range<uint64_t>(), and REDC_non_minimized<uint64_t>() in
    // monty_common.h, for details on what is done here and why it works.
    uint64_t rrax = u_lo;
    uint64_t dummy;
    asm("movq %%rax, %1 \n\t"       /* reg1 = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %1, %%rax \n\t" */ /* rax = u_lo + mn_lo. Sets carry, sum==0 */
          "negq %1 \n\t"            /* Sets carry to u_lo!=0 (same as above) */
        "adcq %[uhi], %%rdx \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        "mov %%rdx, %%rax \n\t"     /* rax = t_hi */
        "subq %[n], %%rdx \n\t"     /* rdx = t_hi - n */
        "cmovaeq %%rdx, %%rax \n\t" /* rax = (t_hi >= n) ? rdx : rax */
        : "+&a"(rrax), "=&r"(dummy)
        : [uhi]"r"(u_hi), [n]"rm"(n), [inv]"rm"(neg_inv_n)
        : "rdx", "cc");
    uint64_t result = rrax;

    HPBC_POSTCONDITION2(result < n);
    HPBC_POSTCONDITION2(result ==
                              montmul_half_range<uint64_t>(x, y, n, neg_inv_n));
    return result;
}
#endif


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
