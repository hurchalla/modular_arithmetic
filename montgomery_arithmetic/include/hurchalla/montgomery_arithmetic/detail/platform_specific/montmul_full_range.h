// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_FULL_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


// The name "full_range" signifies that there are no extra preconditions on
// the value of modulus 'n'.  Although montgomery multiplication always requires
// that modulus 'n' is odd, this function will work for any odd 'n' from the
// full range of possible uint64_t values.


template <typename T>
HURCHALLA_FORCE_INLINE T montmul_full_range(T x, T y, T n, T neg_inv_n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);

    // x<n with y<n  will satisfy montmul_non_minimized's precondition
    // requirement that x*y < n*R.
    bool ovf;
    T prod = montmul_non_minimized(ovf, x, y, n, neg_inv_n);

    // montmul_non_minimized()'s postconditions guarantee the following
    T minimized_result = (ovf || prod >= n) ? static_cast<T>(prod-n) : prod;
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
// This function is an asm version of the template function montmul_full_range()
HURCHALLA_FORCE_INLINE uint64_t montmul_full_range(uint64_t x, uint64_t y,
                                                 uint64_t n, uint64_t neg_inv_n)
{
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);

    // Performing this first multiply outside the asm block lets the compiler
    // potentially optimize it.  E.g. if it can see x == y, or x == 1 or y == 1.
    __uint128_t u = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(y);
    uint64_t u_hi = static_cast<uint64_t>(u >> 64);   // rdx
    uint64_t u_lo = static_cast<uint64_t>(u);         // rax

    // montgomery REDC portion.
    // See montmul_full_range<uint64_t>(), and REDC_non_minimized<uint64_t>() in
    // monty_common.h, for details on what is done here and why it works.
    uint64_t result, reg1;
    asm("movq %%rax, %1 \n\t"       /* reg1 = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %1, %%rax \n\t" */ /* rax = u_lo + mn_lo. Sets carry, rax==0 */
          "xorl %%eax, %%eax \n\t"  /* rax = 0.  CPU will zero upper half too */
          "negq %1 \n\t"            /* u_lo = -u_lo. Sets carry = (u_lo != 0) */
        "movl %%eax, %k1 \n\t"      /* reg1 = 0.  The CPU will zero the upper half too.  For the k modifier see https://gcc.gnu.org/onlinedocs/gcc/Extended-Asm.html#x86Operandmodifiers */
        "adcq %[uhi], %%rdx \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        "cmovaeq %[n], %1 \n\t"     /* reg1 = (t_hi >= u_hi) ? n : 0 */
        "subq %[n], %%rdx \n\t"     /* rdx = t_hi - n */
        "cmovbq %1, %%rax \n\t"     /* rax = (t_hi < n) ? reg1 : 0 */
        "addq %%rdx, %%rax \n\t"    /* result = rax + rdx */
        : "=&a"(result), "=&r"(reg1)
        : "0"(u_lo), [uhi]"r"(u_hi), [n]"r"(n), [inv]"r"(neg_inv_n)
        : "rdx", "cc");

    HPBC_POSTCONDITION2(result < n);
    return result;
}
#endif


}} // end namespace

#endif
