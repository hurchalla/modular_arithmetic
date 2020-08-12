// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTMUL_FULL_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


// The name "fullrange" signifies that there are no extra preconditions on
// the value of modulus 'n'.  Although montgomery multiplication always requires
// that modulus 'n' is odd, this function will work for any odd 'n' from the
// full range of possible type T values.

// The public functions are at the bottom of this file.

namespace detail_fullrange {
// -----------------
// Private Functions
// -----------------


// Primary template.  AV stands for AssemblyVariant (ie. InplaceLowlatencyTag,
// OutofplaceLowlatencyTag, InplaceLowuopsTag, OutofplaceLowuopsTag)
template <typename T, class AV>
struct MontFunctionsFullRange;


// default non-assembly implementation used by montmul_fullrange().
template <typename T>
HURCHALLA_FORCE_INLINE T default_montmul_fullrange(T x, T y, T n, T neg_inv_n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);

    // x<n with y<n  satisfies montmul_non_minimized's precondition
    // requirement that x*y < n*R.
    bool ovf;
    T prod = montmul_non_minimized(ovf, x, y, n, neg_inv_n);

    // montmul_non_minimized()'s postconditions guarantee the following
    T minimized_result = (ovf || prod >= n) ? static_cast<T>(prod-n) : prod;
    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
}

template <typename T>
HURCHALLA_FORCE_INLINE T default_montfmsub_fullrange(T x, T y, T z, T n,
                                                                    T neg_inv_n)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);
    HPBC_PRECONDITION2(z < n);

    // x<n with y<n  satisfies montmul_non_minimized's precondition
    // requirement that x*y < n*R.
    // z<n satisfies montfmsub_non_minimized's precondition that z<n.
    bool ovf;
    T result = montfmsub_non_minimized(ovf, x, y, z, n, neg_inv_n);

    // montfmsub_non_minimized()'s postconditions guarantee the following
    T minimized_result = (ovf || result>=n) ? static_cast<T>(result-n) : result;
    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
}



#if defined(HURCHALLA_ALLOW_INLINE_ASM_MONTMUL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

// InplaceLowlatencyTag version
// This should have: cycles latency 10, fused uops 12
HURCHALLA_FORCE_INLINE std::uint64_t montREDC_fullrange(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                                           InplaceLowlatencyTag)
{
    using std::uint64_t;
    // We require u = (u_hi*R + u_lo) < n*R.  The following guarantees it:
    HPBC_PRECONDITION2(u_hi < n);
    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.

    // montgomery REDC algorithm:
    // See REDC_non_minimized<uint64_t>() in monty_common.h, for details on what
    // is done here and why it works.  We fully minimize the result here.
    uint64_t rrax = u_lo;
    uint64_t rrdx, dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, rax==0 */
          "xorl %%eax, %%eax \n\t"  /* rax = 0.  CPU will zero upper half too */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %[uhi], %%rdx \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        /* "movl %%eax, %k[tmp] \n\t" */  /* tmp = 0.  The CPU will zero the upper half too.  For the k modifier see https://gcc.gnu.org/onlinedocs/gcc/Extended-Asm.html#x86Operandmodifiers */
        "cmovaeq %[n], %%rax \n\t"  /* rax = (t_hi >= u_hi) ? n : 0 */
          "xorl %k[tmp], %k[tmp] \n\t"  /* tmp = 0.  The CPU will zero the upper half too.  For the k modifier see https://gcc.gnu.org/onlinedocs/gcc/Extended-Asm.html#x86Operandmodifiers */
        "subq %[n], %%rdx \n\t"     /* rdx = t_hi - n */
        "cmovaeq %[tmp], %%rax \n\t"  /* rax = (t_hi >= n) ? 0 : rax */
        : "+&a"(rrax), "=&d"(rrdx), [tmp]"=&r"(dummy)
        : [uhi]"r"(u_hi), [n]"r"(n), [inv]"rm"(neg_inv_n)
        : "cc");
    uint64_t result = rrax + rrdx;   // let compiler choose between add/lea

    HPBC_POSTCONDITION2(result < n);
    return result;
}

// OutofplaceLowlatencyTag version
// This should have: cycles latency 10, fused uops 11
HURCHALLA_FORCE_INLINE std::uint64_t montREDC_fullrange(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                                        OutofplaceLowlatencyTag)
{
    using std::uint64_t;
    // We require u = (u_hi*R + u_lo) < n*R.  The following guarantees it:
    HPBC_PRECONDITION2(u_hi < n);
    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.

    // montgomery REDC algorithm:
    // See REDC_non_minimized<uint64_t>() in monty_common.h, for details on what
    // is done here and why it works.  We fully minimize the result here.
    uint64_t rrax = u_lo;
    uint64_t reg = u_hi;
    uint64_t rrdx, dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, rax==0 */
          "xorl %%eax, %%eax \n\t"  /* rax = 0.  CPU will zero upper half too */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %%rdx, %[reg] \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        "cmovbq %[n], %%rax \n\t"   /* subtrahend = (t_hi < u_hi) ? n : 0 */
        "cmpq %[n], %[reg] \n\t"    /* compare (t_hi - n) */
        "cmovaeq %[n], %%rax \n\t"  /* subtrahend = (t_hi>=n) ? n : subtrahnd */
        "subq %%rax, %[reg] \n\t"   /* result = t_hi - subtrahend */
        : [reg]"+&r"(reg), "+&a"(rrax), "=&d"(rrdx), [tmp]"=&r"(dummy)
        : [n]"r"(n), [inv]"rm"(neg_inv_n)
        : "cc");
    uint64_t result = reg;

    HPBC_POSTCONDITION2(result < n);
    return result;
}

// InplaceLowuopsTag version
// This should have: cycles latency 11, fused uops 11
HURCHALLA_FORCE_INLINE std::uint64_t montREDC_fullrange(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                                              InplaceLowuopsTag)
{
    using std::uint64_t;
    // We require u = (u_hi*R + u_lo) < n*R.  The following guarantees it:
    HPBC_PRECONDITION2(u_hi < n);
    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.

    // montgomery REDC algorithm:
    // See REDC_non_minimized<uint64_t>() in monty_common.h, for details on what
    // is done here and why it works.  We fully minimize the result here.
    uint64_t rrax = u_lo;
    uint64_t rrdx, dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
          "addq %[tmp], %%rax \n\t" /* rax = u_lo + mn_lo. Sets carry, rax==0 */
        "adcq %[uhi], %%rdx \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        /* "movl %%eax, %k[tmp] \n\t" */  /* tmp = 0.  The CPU will zero the upper half too.  For the k modifier see https://gcc.gnu.org/onlinedocs/gcc/Extended-Asm.html#x86Operandmodifiers */
        "cmovaeq %[n], %%rax \n\t"  /* rax = (t_hi >= u_hi) ? n : 0 */
          "xorl %k[tmp], %k[tmp] \n\t"  /* tmp = 0.  The CPU will zero the upper half too.  For the k modifier see https://gcc.gnu.org/onlinedocs/gcc/Extended-Asm.html#x86Operandmodifiers */
        "subq %[n], %%rdx \n\t"     /* rdx = t_hi - n */
        "cmovaeq %[tmp], %%rax \n\t"  /* rax = (t_hi >= n) ? 0 : rax */
        : "+&a"(rrax), "=&d"(rrdx), [tmp]"=&r"(dummy)
        : [uhi]"r"(u_hi), [n]"r"(n), [inv]"rm"(neg_inv_n)
        : "cc");
    uint64_t result = rrax + rrdx;   // let compiler choose between add/lea

    HPBC_POSTCONDITION2(result < n);
    return result;
}

// OutofplaceLowuopsTag version
// This should have: cycles latency 11, fused uops 10
HURCHALLA_FORCE_INLINE std::uint64_t montREDC_fullrange(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                                           OutofplaceLowuopsTag)
{
    using std::uint64_t;
    // We require u = (u_hi*R + u_lo) < n*R.  The following guarantees it:
    HPBC_PRECONDITION2(u_hi < n);
    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.

    // montgomery REDC algorithm:
    // See REDC_non_minimized<uint64_t>() in monty_common.h, for details on what
    // is done here and why it works.  We fully minimize the result here.
    uint64_t rrax = u_lo;
    uint64_t reg = u_hi;
    uint64_t rrdx, dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
          "addq %[tmp], %%rax \n\t" /* rax = u_lo + mn_lo. Sets carry, rax==0 */
        "adcq %%rdx, %[reg] \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        "cmovbq %[n], %%rax \n\t"   /* subtrahend = (t_hi < u_hi) ? n : 0 */
        "cmpq %[n], %[reg] \n\t"    /* compare (t_hi - n) */
        "cmovaeq %[n], %%rax \n\t"  /* subtrahend = (t_hi>=n) ? n : subtrahnd */
        "subq %%rax, %[reg] \n\t"   /* result = t_hi - subtrahend */
        : [reg]"+&r"(reg), "+&a"(rrax), "=&d"(rrdx), [tmp]"=&r"(dummy)
        : [n]"r"(n), [inv]"rm"(neg_inv_n)
        : "cc");
    uint64_t result = reg;

    HPBC_POSTCONDITION2(result < n);
    return result;
}


template <class AV>
struct MontFunctionsFullRange<std::uint64_t, AV>
{
  using T = std::uint64_t;

  // Expected latency and uops for mul() (includes REDC latency/uops), for AV ==
  //   InplaceLowlatencyTag     cycles latency 13, fused uops 15
  //   OutofplaceLowlatencyTag  cycles latency 13, fused uops 14
  //   InplaceLowuopsTag        cycles latency 14, fused uops 14
  //   OutofplaceLowuopsTag     cycles latency 14, fused uops 13
  static HURCHALLA_FORCE_INLINE T mul(T x, T y, T n, T neg_inv_n)
  {
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);
    // Performing this multiply outside of an asm block lets the compiler
    // potentially optimize it.  E.g. if it can see x == y, or x == 1 or y == 1.
    __uint128_t u = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(y);
    T u_hi = static_cast<T>(u >> 64);   // rdx
    T u_lo = static_cast<T>(u);         // rax

    T result = montREDC_fullrange(u_hi, u_lo, n, neg_inv_n, AV());

    HPBC_POSTCONDITION2(result < n);
    HPBC_POSTCONDITION2(result== default_montmul_fullrange(x, y, n, neg_inv_n));
    return result;
  }

  // TODO need to refer to proof that the integrated modular subtract works
  // Expected latency and uops for fmsub() (includes REDC latency/uops), AV ==
  //   InplaceLowlatencyTag     cycles latency 13, fused uops 17
  //   OutofplaceLowlatencyTag  cycles latency 13, fused uops 16
  //   InplaceLowuopsTag        cycles latency 14, fused uops 16
  //   OutofplaceLowuopsTag     cycles latency 14, fused uops 15
  static HURCHALLA_FORCE_INLINE T fmsub(T x, T y, T z, T n, T neg_inv_n)
  {
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);
    HPBC_PRECONDITION2(z < n);
    // Performing this multiply outside of an asm block lets the compiler
    // potentially optimize it.  E.g. if it can see x == y, or x == 1 or y == 1.
    __uint128_t u = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(y);
    T u_hi = static_cast<T>(u >> 64);   // rdx
    T u_lo = static_cast<T>(u);         // rax

    T rrdx = u_hi;
    T u_hi2;
    asm("subq %[z], %%rdx \n\t"             /* rdx = u_hi - z */
        "leaq (%%rdx, %[n]), %[uhi2] \n\t"  /* uhi2 = u_hi - z + n */
        "cmovaeq %%rdx, %[uhi2] \n\t"       /* uhi2 = (u_hi>=z) ? rdx : uhi2 */
        : "+&d"(rrdx), [uhi2]"=&r"(u_hi2)
        : [n]"r"(n), [z]"rm"(z)
        : "cc");
    T result = montREDC_fullrange(u_hi2, u_lo, n, neg_inv_n, AV());

    HPBC_POSTCONDITION2(result < n);
    HPBC_POSTCONDITION2(result == default_montfmsub_fullrange(x, y, z, n,
                                                                    neg_inv_n));
    return result;
  }
};

#endif   // defined(HURCHALLA_ALLOW_INLINE_ASM_MONTMUL) &&
         // defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// -------------------------  END OF ASM ------------------------------------


// The default (non-asm) implementations for montmul and montfmsub full range.
template <typename T, class AV>
struct MontFunctionsFullRange
{
    static HURCHALLA_FORCE_INLINE T mul(T x, T y, T n, T neg_inv_n)
    {
        return default_montmul_fullrange(x, y, n, neg_inv_n);
    }
    static HURCHALLA_FORCE_INLINE T fmsub(T x, T y, T z, T n, T neg_inv_n)
    {
        return default_montfmsub_fullrange(x, y, z, n, neg_inv_n);
    }
};


}  // end namespace detail_fullrange (end of the private classes/functions)

// -----------------
// Public Functions
// -----------------


template <typename T, class AV=OutofplaceLowlatencyTag>
HURCHALLA_FORCE_INLINE T montmul_fullrange(T x, T y, T n, T neg_inv_n)
{
    // Delegate to a static class function, since it's impossible to partially
    // specialize a function (instead, we partially specialize the class).
    T result = detail_fullrange::MontFunctionsFullRange<T, AV>::mul(x, y, n,
                                                                     neg_inv_n);
    HPBC_POSTCONDITION2(result < n);
    return result;
}

template <typename T, class AV=OutofplaceLowlatencyTag>
HURCHALLA_FORCE_INLINE T montfmsub_fullrange(T x, T y, T z, T n, T neg_inv_n)
{
    // Delegate to a static class function, since it's impossible to partially
    // specialize a function (instead, we partially specialize the class).
    T result = detail_fullrange::MontFunctionsFullRange<T, AV>::fmsub(x, y, z,
                                                                  n, neg_inv_n);
    HPBC_POSTCONDITION2(result < n);
    return result;
}


}} // end namespace

#endif
