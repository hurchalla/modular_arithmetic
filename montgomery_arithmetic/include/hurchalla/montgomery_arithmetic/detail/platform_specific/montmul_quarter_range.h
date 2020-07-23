// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_QUARTER_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_QUARTER_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#  pragma warning(disable : 4189)
#endif

namespace hurchalla { namespace montgomery_arithmetic {


// The name "quarter_range" signifies that modulus 'n' must be less than R/4,
// where  R = 2^(ma_numeric_limits<T>::digits).  For example, if T is uint64_t
// then R = 2^64 and R/4 == 2^62, and thus we would require  n < 2^62.


// montmul_quarter_range() requires/allows an unusual input range:  we allow
// 0 <= x < 2*n,  and  0 <= y < 2*n.
// Similarly, the return value range will be  0 <= returnValue < 2*n.
// Obviously neither inputs nor outputs necessarily belong to the minimal
// residue class modulo n -- i.e. they might not be fully reduced, modulo n.
// The preconditions/postconditions in impl_montmul_non_minimized()  (from
// monty_common.h) have accompanying proofs that this works, but for more
// details, see also section 5 of the paper "Montgomery's Multiplication
// Technique: How to Make It Smaller and Faster"
// https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps
template <typename T>
HURCHALLA_FORCE_INLINE T montmul_quarter_range(T x, T y, T n, T neg_inv_n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(n < Rdiv4);
    }
    HPBC_PRECONDITION2(x < 2*n);
    HPBC_PRECONDITION2(y < 2*n);

    // Since x<2*n and y<2*n, we know x*y < 4*n*n, and since the modulus
    // n < R/4, we know  x*y < 4*n*R/4 == n*R.
    // This satisfies montmul_non_minimized's precondition of x*y < n*R.
    bool ovf;
    T prod = montmul_non_minimized(ovf, x, y, n, neg_inv_n);

    // Since we know that modulus n < R/4, the postconditions of
    // montmul_non_minimized() guarantee  prod < 2*n.
    HPBC_POSTCONDITION2(prod < 2*n);
    return prod;
}


// -------- PLATFORM SPECIFIC nontemplate overloads ----------

// Note: a nontemplate function overload gets first priority for being called
// (see http://www.gotw.ca/publications/mill17.htm ), when both the nontemplate
// function and the generic template function match the caller's provided
// argument type(s).


#if defined(HURCHALLA_ALLOW_INLINE_ASM_MONTMUL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// This function is an asm version of template function montmul_quarter_range().
HURCHALLA_FORCE_INLINE uint64_t montmul_quarter_range(uint64_t x, uint64_t y,
                                                 uint64_t n, uint64_t neg_inv_n)
{
    HPBC_PRECONDITION2(n < (static_cast<uint64_t>(1) << 62));
    HPBC_PRECONDITION2(x < 2*n);
    HPBC_PRECONDITION2(y < 2*n);

    // Performing this first multiply outside the asm block lets the compiler
    // potentially optimize it.  E.g. if it can see x == y, or x == 1 or y == 1.
    __uint128_t u = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(y);
    uint64_t u_hi = static_cast<uint64_t>(u >> 64);   // rdx
    uint64_t u_lo = static_cast<uint64_t>(u);         // rax

    // montgomery REDC portion.
    // See montmul_quarter_range<uint64_t>(), and REDC_non_minimized<uint64_t>() in
    // monty_common.h, for details on what is done here and why it works.
    uint64_t result, reg1;
    asm("movq %%rax, %1 \n\t"       /* reg1 = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %1, %%rax \n\t" */ /* rax = u_lo + mn_lo. Sets carry, sum==0 */
          "movq %[uhi], %%rax \n\t"
          "negq %1 \n\t"            /* Sets carry to u_lo!=0. Carry eqls above*/
        "adcq %%rdx, %%rax \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        : "=&a"(result), "=&r"(reg1)
        : "0"(u_lo), [uhi]"r"(u_hi), [n]"r"(n), [inv]"r"(neg_inv_n)
        : "rdx", "cc");

    HPBC_POSTCONDITION2(result < 2*n);
    HPBC_POSTCONDITION2(result ==
                           montmul_quarter_range<uint64_t>(x, y, n, neg_inv_n));
    return result;
}
#endif


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
