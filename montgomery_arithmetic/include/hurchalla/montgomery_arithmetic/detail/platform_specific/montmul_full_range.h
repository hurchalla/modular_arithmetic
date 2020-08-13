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


// The public functions are at the bottom of this file.

namespace mf_detail {
// -----------------
// Private Functions
// -----------------


// MT stands for the Montgomery type (ie. FullrangeTag, HalfrangeTag,
// QuarterrangeTag)
// Primary template.  AV stands for AssemblyVariant (ie. InplaceLowlatencyTag,
// OutofplaceLowlatencyTag, InplaceLowuopsTag, OutofplaceLowuopsTag)
template <typename T, class MT, class AV>
struct MontFunctions;



// The name "fullrange" signifies that there are no extra preconditions on
// the value of modulus 'n'.  Although montgomery multiplication always requires
// that modulus 'n' is odd, this function will work for any odd 'n' from the
// full range of possible type T values.
template <typename T> HURCHALLA_FORCE_INLINE
void check_montmul_preconditions(T x, T y, T n, FullrangeTag)
{
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(x < n);
    HPBC_PRECONDITION2(y < n);
    // Prior to performing REDC any multiply sets  u = x*y.
    // Thus x<n with y<n  satisfies the montgomery REDC algorithm's requirement 
    // that the REDC input u < n*R.
}

// The name "halfrange" signifies that modulus 'n' must be less than R/2, where
// R = 2^(ma_numeric_limits<T>::digits).  For example, if T is uint64_t then
// R = 2^64 and R/2 == 2^63, and thus we would require  n < 2^63.
template <typename T> HURCHALLA_FORCE_INLINE
void check_montmul_preconditions(T x, T y, T n, HalfrangeTag)
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
    // Prior to performing REDC any multiply sets  u = x*y.
    // Thus x<n with y<n  satisfies the montgomery REDC algorithm's requirement 
    // that the REDC input u < n*R.
}

// The name "quarter_range" signifies that modulus 'n' must be less than R/4,
// where  R = 2^(ma_numeric_limits<T>::digits).  For example, if T is uint64_t
// then R = 2^64 and R/4 == 2^62, and thus we would require  n < 2^62.
//
// montmul_quarter_range() requires/allows an unusual input range:  we allow
// 0 <= x < 2*n,  and  0 <= y < 2*n.
// Similarly, the return value range will be  0 <= returnValue < 2*n.
// Obviously neither inputs nor outputs necessarily belong to the minimal
// residue class modulo n -- i.e. they might not be fully reduced, modulo n.
// The preconditions/postconditions in montmul_non_minimized()  (from
// monty_common.h) have accompanying proofs that this works, but for more
// details, see also section 5 of the paper "Montgomery's Multiplication
// Technique: How to Make It Smaller and Faster"
// https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps
template <typename T> HURCHALLA_FORCE_INLINE
void check_montmul_preconditions(T x, T y, T n, QuarterrangeTag)
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
    // Since x<2*n and y<2*n, we know x*y < 4*n*n, and since one of our
    // preconditions was n < R/4, we know  x*y < 4*n*R/4 == n*R.
    // Prior to performing REDC any multiply sets  u = x*y.  Since x*y < n*R, we
    // have u = x*y < n*R, which satisfies the montgomery REDC algorithm's 
    // requirement that the REDC input u < n*R.
}



// default non-assembly implementations
template <typename T> HURCHALLA_FORCE_INLINE
T default_montmul(T x, T y, T n, T neg_inv_n, FullrangeTag)
{
    check_montmul_preconditions(x, y, n, FullrangeTag());
    bool ovf;
    T product = montmul_non_minimized(ovf, x, y, n, neg_inv_n);
    // montmul_non_minimized()'s postconditions guarantee the following
    T minimized_prod= (ovf || product>=n) ? static_cast<T>(product-n) : product;
    HPBC_POSTCONDITION2(minimized_prod < n);
    return minimized_prod;
}
template <typename T> HURCHALLA_FORCE_INLINE
T default_montfmsub(T x, T y, T z, T n, T neg_inv_n, FullrangeTag)
{
    check_montmul_preconditions(x, y, n, FullrangeTag());
    HPBC_PRECONDITION2(z < n);   // z must always be canonical.
    bool ovf;
    T product = montfmsub_non_minimized(ovf, x, y, z, n, neg_inv_n);
    // montfmsub_non_minimized()'s postconditions guarantee the following
    T minimized_prod= (ovf || product>=n) ? static_cast<T>(product-n) : product;
    HPBC_POSTCONDITION2(minimized_prod < n);
    return minimized_prod;
}

template <typename T> HURCHALLA_FORCE_INLINE
T default_montmul(T x, T y, T n, T neg_inv_n, HalfrangeTag)
{
    check_montmul_preconditions(x, y, n, HalfrangeTag());
    bool ovf;
    T product = montmul_non_minimized(ovf, x, y, n, neg_inv_n);

    // Since the preconditions for HalfrangeTag require  n < R/2,  we know by
    // montmul_non_minimized()'s postconditions that ovf must == false.
    HPBC_ASSERT2(ovf == false);
    // Since ovf == false, montmul_non_minimized()'s postconditions guarantee
    T minimized_product = (product >= n) ? static_cast<T>(product-n) : product;
    HPBC_POSTCONDITION2(minimized_product < n);
    return minimized_product;
}
template <typename T> HURCHALLA_FORCE_INLINE
T default_montfmsub(T x, T y, T z, T n, T neg_inv_n, HalfrangeTag)
{
    check_montmul_preconditions(x, y, n, HalfrangeTag());
    HPBC_PRECONDITION2(z < n);   // z must always be canonical.
    bool ovf;
    T product = montfmsub_non_minimized(ovf, x, y, z, n, neg_inv_n);

    // Since the preconditions for HalfrangeTag require  n < R/2,  we know by
    // montfmsub_non_minimized()'s postconditions that ovf must == false.
    HPBC_ASSERT2(ovf == false);
    // Since ovf == false, montfmsub_non_minimized()'s postconditions guarantee
    T minimized_product = (product >= n) ? static_cast<T>(product-n) : product;
    HPBC_POSTCONDITION2(minimized_product < n);
    return minimized_product;
}

template <typename T> HURCHALLA_FORCE_INLINE
T default_montmul(T x, T y, T n, T neg_inv_n, QuarterrangeTag)
{
    check_montmul_preconditions(x, y, n, QuarterrangeTag());
    bool ovf;
    T product = montmul_non_minimized(ovf, x, y, n, neg_inv_n);

    // Since the preconditions for QuarterrangeTag require  n < R/4,  we know by
    // montmul_non_minimized()'s postconditions that  product < 2*n.
    // QuarterrangeTag uses non-minimized inputs/outputs.  0<=product<2*n is ok.
    HPBC_POSTCONDITION2(product < 2*n);
    return product;
}
template <typename T> HURCHALLA_FORCE_INLINE
T default_montfmsub(T x, T y, T z, T n, T neg_inv_n, QuarterrangeTag)
{
    check_montmul_preconditions(x, y, n, QuarterrangeTag());
    HPBC_PRECONDITION2(z < n);   // z must always be canonical.
    bool ovf;
    T product = montfmsub_non_minimized(ovf, x, y, z, n, neg_inv_n);

    // Since the preconditions for QuarterrangeTag require  n < R/4,  we know by
    // montfmsub_non_minimized()'s postconditions that  product < 2*n.
    // QuarterrangeTag uses non-minimized inputs/outputs.  0<=product<2*n is ok.
    HPBC_POSTCONDITION2(product < 2*n);
    return product;
}



#if defined(HURCHALLA_ALLOW_INLINE_ASM_MONTMUL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

// FullrangeTag InplaceLowlatencyTag version
// This should have: cycles latency 10, fused uops 12
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                             FullrangeTag, InplaceLowlatencyTag)
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

// FullrangeTag OutofplaceLowlatencyTag version
// This should have: cycles latency 10, fused uops 11
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                          FullrangeTag, OutofplaceLowlatencyTag)
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

// FullrangeTag InplaceLowuopsTag version
// This should have: cycles latency 11, fused uops 11
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                                FullrangeTag, InplaceLowuopsTag)
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

// FullrangeTag OutofplaceLowuopsTag version
// This should have: cycles latency 11, fused uops 10
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                             FullrangeTag, OutofplaceLowuopsTag)
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



// HalfrangeTag PrivateInplaceTag version (covers InplaceLowlatencyTag and
// InplaceLowuopsTag)
// This should have: cycles latency 9, fused uops 9
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                                HalfrangeTag, PrivateInplaceTag)
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
    uint64_t dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, sum==0 */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %[uhi], %%rdx \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        "mov %%rdx, %%rax \n\t"     /* rax = t_hi */
        "subq %[n], %%rdx \n\t"     /* rdx = t_hi - n */
        "cmovaeq %%rdx, %%rax \n\t" /* rax = (t_hi >= n) ? rdx : rax */
        : "+&a"(rrax), [tmp]"=&r"(dummy)
        : [uhi]"r"(u_hi), [n]"rm"(n), [inv]"rm"(neg_inv_n)
        : "rdx", "cc");
    uint64_t result = rrax;

    HPBC_POSTCONDITION2(result < n);
    return result;
}

// HalfrangeTag PrivateOutofplaceTag version (covers OutofplaceLowlatencyTag and
// OutofplaceLowuopsTag)
// This should have: cycles latency 9, fused uops 9
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                             HalfrangeTag, PrivateOutofplaceTag)
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
    uint64_t dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, sum==0 */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %%rdx, %[reg] \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        "mov %[reg], %%rax \n\t"    /* rax = t_hi */
        "subq %[n], %[reg] \n\t"    /* reg = t_hi - n */
        "cmovbq %%rax, %[reg] \n\t" /* reg = (t_hi < n) ? rax : reg */
        : [reg]"+&r"(reg), "+&a"(rrax), [tmp]"=&r"(dummy)
        : [n]"rm"(n), [inv]"rm"(neg_inv_n)
        : "rdx", "cc");
    uint64_t result = reg;

    HPBC_POSTCONDITION2(result < n);
    return result;
}



// QuarterrangeTag PrivateInplaceTag version (covers InplaceLowlatencyTag and
// InplaceLowuopsTag)
// This should have: cycles latency 7, fused uops 7
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                             QuarterrangeTag, PrivateInplaceTag)
{
    using std::uint64_t;
    // We require u = (u_hi*R + u_lo) < n*R.  The following guarantees it:
    HPBC_PRECONDITION2(u_hi < n);
    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.

    // montgomery REDC algorithm:
    // See REDC_non_minimized<uint64_t>() in monty_common.h, for details on what
    // is done here and why it works.  We do *not* minimize the result for
    // QuarterrangeTag.
    uint64_t rrax = u_lo;
    uint64_t dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, sum==0 */
          "movq %[uhi], %%rax \n\t" /* rax = u_hi */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %%rdx, %%rax \n\t"    /* t_hi = addcarry(u_hi, mn_hi) */
        : "+&a"(rrax), [tmp]"=&r"(dummy)
        : [uhi]"r"(u_hi), [n]"rm"(n), [inv]"rm"(neg_inv_n)
        : "rdx", "cc");
    uint64_t result = rrax;

    HPBC_POSTCONDITION2(result < 2*n);
    return result;
}

// QuarterrangeTag PrivateOutofplaceTag version (covers OutofplaceLowlatencyTag
// and OutofplaceLowuopsTag)
// This should have: cycles latency 7, fused uops 6
HURCHALLA_FORCE_INLINE std::uint64_t montREDC(std::uint64_t u_hi,
                   std::uint64_t u_lo, std::uint64_t n, std::uint64_t neg_inv_n,
                                          QuarterrangeTag, PrivateOutofplaceTag)
{
    using std::uint64_t;
    // We require u = (u_hi*R + u_lo) < n*R.  The following guarantees it:
    HPBC_PRECONDITION2(u_hi < n);
    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.

    // montgomery REDC algorithm:
    // See REDC_non_minimized<uint64_t>() in monty_common.h, for details on what
    // is done here and why it works.  We do *not* minimize the result for
    // QuarterrangeTag.
    uint64_t rrax = u_lo;
    uint64_t reg = u_hi;
    uint64_t dummy;
    asm("movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, sum==0 */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %%rdx, %[reg] \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        : [reg]"+&r"(reg), "+&a"(rrax), [tmp]"=&r"(dummy)
        : [n]"rm"(n), [inv]"rm"(neg_inv_n)
        : "rdx", "cc");
    uint64_t result = reg;

    HPBC_POSTCONDITION2(result < 2*n);
    return result;
}




// MT stands for the Montgomery type (ie. FullrangeTag, HalfrangeTag,
// QuarterrangeTag)
// AV stands for AssemblyVariant (ie. InplaceLowlatencyTag,
// OutofplaceLowlatencyTag, InplaceLowuopsTag, OutofplaceLowuopsTag)

template <class MT, class AV>
struct MontFunctions<std::uint64_t, MT, AV>
{
  using T = std::uint64_t;

  // TODO stats for other MT tags

  // Expected latency and uops for mul() (includes REDC latency/uops), for
  // MT == FullrangeTag and AV ==
  //   InplaceLowlatencyTag     cycles latency 13, fused uops 15
  //   OutofplaceLowlatencyTag  cycles latency 13, fused uops 14
  //   InplaceLowuopsTag        cycles latency 14, fused uops 14
  //   OutofplaceLowuopsTag     cycles latency 14, fused uops 13
  static HURCHALLA_FORCE_INLINE T mul(T x, T y, T n, T neg_inv_n)
  {
    check_montmul_preconditions(x, y, n, MT());
    // Performing this multiply outside of an asm block lets the compiler
    // potentially optimize it.  E.g. if it can see x == y, or x == 1 or y == 1.
    __uint128_t u = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(y);
    T u_hi = static_cast<T>(u >> 64);   // rdx
    T u_lo = static_cast<T>(u);         // rax

    T result = montREDC(u_hi, u_lo, n, neg_inv_n, MT(), AV());

    HPBC_POSTCONDITION2(result == default_montmul(x, y, n, neg_inv_n, MT()));
    return result;
  }

  // TODO need to refer to proof that the integrated modular subtract works
  // TODO stats for other MT tags

  // Expected latency and uops for fmsub() (includes REDC latency/uops), for
  // MT == FullrangeTag and AV ==
  //   InplaceLowlatencyTag     cycles latency 13, fused uops 17
  //   OutofplaceLowlatencyTag  cycles latency 13, fused uops 16
  //   InplaceLowuopsTag        cycles latency 14, fused uops 16
  //   OutofplaceLowuopsTag     cycles latency 14, fused uops 15
  static HURCHALLA_FORCE_INLINE T fmsub(T x, T y, T z, T n, T neg_inv_n)
  {
    check_montmul_preconditions(x, y, n, MT());
    HPBC_PRECONDITION2(z < n);  // z must always be canonical.
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
    T result = montREDC(u_hi2, u_lo, n, neg_inv_n, MT(), AV());

    HPBC_POSTCONDITION2(result == default_montfmsub(x,y,z,n,neg_inv_n,MT()));
    return result;
  }
};

#endif   // defined(HURCHALLA_ALLOW_INLINE_ASM_MONTMUL) &&
         // defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
// -------------------------  END OF ASM ------------------------------------


// The default (non-asm) implementations for montmul and montfmsub full range.
template <typename T, class MT, class AV>
struct MontFunctions
{
    static HURCHALLA_FORCE_INLINE T mul(T x, T y, T n, T neg_inv_n)
    {
        return default_montmul(x, y, n, neg_inv_n, MT());
    }
    static HURCHALLA_FORCE_INLINE T fmsub(T x, T y, T z, T n, T neg_inv_n)
    {
        return default_montfmsub(x, y, z, n, neg_inv_n, MT());
    }
};


}  // end namespace mf_detail (end of the private classes/functions)

// -----------------
// Public Functions
// -----------------


template <typename T, class MT, class AV=OutofplaceLowlatencyTag>
HURCHALLA_FORCE_INLINE T impl_montmul(T x, T y, T n, T neg_inv_n)
{
    mf_detail::check_montmul_preconditions(x, y, n, MT());
    // Delegate to a static class function, since it's impossible to partially
    // specialize a function (instead, we partially specialize the class).
    T result = mf_detail::MontFunctions<T,MT,AV>::mul(x, y, n, neg_inv_n);
    return result;
}

template <typename T, class MT, class AV=OutofplaceLowlatencyTag>
HURCHALLA_FORCE_INLINE T impl_montfmsub(T x, T y, T z, T n, T neg_inv_n)
{
    mf_detail::check_montmul_preconditions(x, y, n, MT());
    // Delegate to a static class function, since it's impossible to partially
    // specialize a function (instead, we partially specialize the class).
    T result = mf_detail::MontFunctions<T,MT,AV>::fmsub(x, y, z, n, neg_inv_n);
    return result;
}


}} // end namespace

#endif


