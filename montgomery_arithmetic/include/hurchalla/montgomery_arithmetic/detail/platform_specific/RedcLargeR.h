// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_REDC_LARGE_R_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_REDC_LARGE_R_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


namespace detail_redc_large {
// -----------------
// Private Functions
// -----------------

// Generic Montgomery REDC algorithm

// This function implements the REDC algorithm as described at
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
// We use the variable name "u" instead of "T" (from the link above) since
// "T" in C++ is conventionally reserved for use as a template parameter name.
// We also use "n" instead of "N", and "neg_inv_n" instead of "N′" (N with a
// prime symbol).  The constant "R" remains the same, and represents the value
// R = 2^(ma::ma_numeric_limits<T>::digits).  As an example, if T is uint64_t,
// then R = 2^64.
// The comments below explain additional changes from the wikipedia article.
//
// This function's name is "REDC_non_minimized" to reflect the fact that the
// value it returns is not minimized to the minimum residual class mod n. (see
// the starred* comment under Postcondition #1 for more info).
template <typename T> HURCHALLA_FORCE_INLINE
T REDC_non_minimized(bool& ovf, T u_hi, T u_lo, T n, T neg_inv_n)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T>::is_modulo, "");
    // We require T is at least as large as unsigned int.
    // If T is a smaller type than that, we'll need a different function that
    // will be both more efficient, and protected from surprises and undefined
    // behavior from the tricky unsigned integral promotion rules in C++.  See
    // https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/
    static_assert(ma::ma_numeric_limits<T>::digits >=
                  ma::ma_numeric_limits<unsigned int>::digits, "");

    // Precondition #1:  We require the precondition  u < n*R.  Or elaborated,
    // u == u_hi*R + u_lo < n*R.
    // If u_hi < n:  then u_hi+1 <= n, and u_hi*R + R <= n*R.  Since u_lo < R,
    //   u == u_hi*R + u_lo < u_hi*R + R <= n*R.  We would have  u < n*R, and
    //   so u_hi < n is always sufficient to satisfy the precondition.
    // If u_hi >= n:  then u_hi*R >= n*R, and u == u_hi*R + u_lo >= n*R, which
    //   fails the precondition.
    // Thus u_hi < n is sufficient and necessary to satisfy the precondition.
    HPBC_PRECONDITION2(u_hi < n);

    // assert(n * neg_inv_n ≡ -1 (mod R))
    HPBC_PRECONDITION2(n * neg_inv_n == static_cast<T>(0) - static_cast<T>(1));
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T m = u_lo * neg_inv_n;  // computes  m = (u * neg_inv_n) % R

    T mn_lo;
    T mn_hi = unsigned_multiply_to_hilo_product(&mn_lo, m, n);

    // mn = m*n.  Since m=(u_lo*neg_inv_n)%R, we know m < R, and thus  mn < R*n.
    // Therefore mn == mn_hi*R + mn_lo < R*n, and mn_hi*R < R*n - mn_lo <= R*n,
    // and thus  mn_hi < n.
        // *** Assertion #1 ***
    HPBC_ASSERT2(mn_hi < n);
    // Since mn_hi < n and n <= R-1,  mn_hi < R-1, or  mn_hi <= R-2.
    // (Note: R-2 == static_cast<T>(0) - static_cast<T>(2))
        // *** Assertion #2 ***  (mn_hi <= R-2)
    HPBC_ASSERT2(mn_hi <= static_cast<T>(0) - static_cast<T>(2));

    // Compute (u + mn)/R :

    T t_hi = u_hi + mn_hi;   // computes  t_hi = (u_hi + mn_hi) % R

    // We do not need to explicitly perform the low part addition (u_lo + mn_lo)
    // because the REDC algorithm guarantees (u_lo + mn_lo) % R == 0.  The only
    // question we have is whether (u_lo + mn_lo) has a carry to the high part.
    HPBC_ASSERT2(u_lo + mn_lo == 0);
    // This low part addition will always carry over into the high part of the
    // sum, unless u_lo and mn_lo are both equal to 0.  Thus, if u_lo != 0, this
    // addition carries.  If u_lo == 0, then since mn_lo < R, we know  mn_lo
    // would need to equal 0 in order to have the low part addition == 0.  Thus,
    // if (u_lo == 0), this addition will not carry.
    // To know if we carry, we simply check (u_lo != 0):
    t_hi += (u_lo != 0);
    // The above computed additions were t_hi = (u_hi + mn_hi + (u_lo != 0)) % R
    // Let  sum = u_hi + mn_hi + (u_lo != 0).  We know (u_lo != 0) <= 1, and
    // assertion #2 shows mn_hi <= R-2, so  sum <= u_hi + R-2 + 1 == u_hi + R-1.
    // Or more simply, sum < u_hi + R.  For some integer k>=0, t_hi + k*R == sum
    // Assume k>=2:  then t_hi + 2*R <= t_hi + k*R == sum < u_hi + R.  More
    // simply, t_hi + R < u_hi.  But since u_hi < R, this means t_hi < 0  which
    // is impossible.  Therefore k<2, and thus k==0 or k==1.
    // If k==1, then  t_hi + R == sum < u_hi + R,  or more simply, t_hi < u_hi.
    // If k==0, then  t_hi == sum == u_hi + mn_hi + (u_lo != 0).  Since
    // mn_hi >= 0 and (u_lo != 0) >= 0,  we would have  t_hi >= u_hi.
    // The contrapositives of these show: if t_hi < u_hi then k!=0, and if
    // t_hi >= u_hi then k!=1.  Since k== 0 or 1,  t_hi < u_hi implies k==1, and
    // t_hi >= u_hi  implies k==0.  This means we can get the value of k using
    // the conditional  k = (t_hi < u_hi).
    // We can view k as an overflow or carry flag for the computed additions in
    // t_hi (we would view the carry as going to a "super high" part of 't').
    // Thus,  bool overflow = k = (t_hi < u_hi)
    ovf = (t_hi < u_hi);

    // The precondition  u < n*R, along with assertion #1 of  mn < n*R,  show
    // that  u + mn < 2*n*R, and thus  t = u + mn < 2*n*R, or  t < 2*n*R.
    // We know t == ovf*R*R + t_hi*R + t_lo, and since we know  t_lo == 0,
    // t == ovf*R*R + t_hi*R.  We need the value q = t/R == ovf*R + t_hi, and
    // since t < 2*n*R, q = t/R < 2*n.  Our desired final result will need to
    // satisfy 0 <= result < n, so we will subtract n from q if q >= n.  Since
    // q < 2*n, we know  q - n < n, showing that a single subtraction at most
    // might be needed to compute the final result.
    // If ovf == 0:
    //   We would know  q = t/R == ovf*R + t_hi == t_hi.  We can't know if we
    //   would have q = t_hi >= n, so we would need to test for it.  If true,
    //   we would perform the subtraction t_hi - n  to get the final result.
    //   Otherwise  t_hi < n  and thus  t_hi  is the final result.
    // If ovf == 1:
    //   We would know  q = t/R == ovf*R + t_hi == R + t_hi.  Since R > n,
    //   we would have  q = R + t_hi > n.  Thus we would need to do a single
    //   subtraction (as proven above) to get the final result.  We would have
    //   result = R + t_hi - n.  Since R + t_hi - n == q - n < n < R  and
    //   since trivially we can see that  R + t_hi - n > 0,  we know 
    //   R + t_hi - n == (R + t_hi - n)%R == (t_hi - n)%R,  which we can compute
    //   with type T variables by  result = t_hi - n  (since type T arithmetic
    //   is implicitly mod R).
    // All this translates into -
    //
    // Postcondition #1 -
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (ovf || t_hi >= n) ? (t_hi - n) : t_hi;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    // * Aside from this postcondition, we do not actually compute the final
    // minimized residual mod n result, because some Montgomery Forms are
    // constrained in ways that allow simpler and more efficient computation of
    // the minimized result.  For example, in some forms, ovf is always false.
    // We allow the caller of this function to compute its final minimized
    // result in the most efficient manner available to it.


    // Postcondition #2 - If u_hi == 0, then ovf == false.
    // ---------------------------------------------------
    // Since ovf = (t_hi < u_hi),  u_hi == 0 would mean ovf = (t_hi < 0), which
    // is always false.
    HPBC_POSTCONDITION2((u_hi == 0) ? ovf == false : true);

    // Postcondition #3 - If u_hi == 0 and u_lo < n, then t_hi < n.
    // ------------------------------------------------------------
    // By definition, m*n == mn_hi*R + mn_lo, and since m is type T,  m <= R-1.
    // Thus mn_hi*R + mn_lo == m*n <= (R-1)*n.  Adding u_lo to both sides we get
    // mn_hi*R + (mn_lo + u_lo) <= (R-1)*n + u_lo.  We know from comments above
    // that if (u_lo == 0), then (u_lo + mn_lo) == 0, and if (u_lo != 0), then
    // (u_lo + mn_lo) == R.  Thus, (u_lo + mn_lo) == (u_lo != 0)*R.  And so
    // mn_hi*R + (u_lo != 0)*R == mn_hi*R + (mn_lo + u_lo) <= (R-1)*n + u_lo.
    // This note specifies  u_lo < n, and so we know  u_lo <= n-1, and thus,
    // mn_hi*R + (u_lo != 0)*R <= (R-1)*n + n-1 == R*n - 1.  Therefore,
    // mn_hi*R + (u_lo != 0)*R < R*n.  Dividing by R,  mn_hi + (u_lo != 0) < n.
    // We know that  t_hi = (u_hi + mn_hi + (u_lo != 0)) % R.  Since this case
    // specifies u_hi == 0, we know  t_hi = (mn_hi + (u_lo != 0)) % R.  And
    // since we know  0 <= mn_hi + (u_lo != 0) < n < R,  we would have
    // t_hi == mn_hi + (u_lo != 0),  and thus,  t_hi < n.
    HPBC_POSTCONDITION2((u_hi == 0 && u_lo < n) ? t_hi < n : true);

    // Postcondition #4 - If n < R/2, then ovf == false  and  t_hi < 2*n.
    // ------------------------------------------------------------------
    // By definition, m*n == mn_hi*R + mn_lo.  This case specifies n < R/2, and
    // since m is type T, we know m < R.  Thus, mn_hi*R + mn_lo == m*n < R*R/2.
    // We know  mn_hi*R <= mn_hi*R + mn_lo < R*R/2, and thus  mn_hi < R/2,  and
    // by extension  mn_hi <= R/2 - 1.  We know (u_lo != 0) <= 1, and by this
    // function's precondition of (u_hi < n) we know  u_hi < n <= R/2 - 1.  Thus
    // u_hi + mn_hi + (u_lo != 0) <= (R/2 - 1) + (R/2 - 1) + 1 == R-1.  And thus
    // 0 <= u_hi + mn_hi + (u_lo != 0) <= R-1 < R.  This lets us see that
    // t_hi = (u_hi + mn_hi + (u_lo != 0)) % R == u_hi + mn_hi + (u_lo != 0).
    // We know t_hi == u_hi + mn_hi + (u_lo != 0) >= u_hi, and so  t_hi >= u_hi.
    // Since  ovf = (t_hi < u_hi),  ovf must be false.
    // From the comments preceding Postcondition #1, we know q == ovf*R + t_hi
    // and q < 2*n.  Since we have ovf == false, q == t_hi, and thus t_hi < 2*n.
    // [ These comments were inspired by ideas in section 5 of "Montgomery's
    // Multiplication Technique: How to Make It Smaller and Faster", at
    // https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps ]
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (ovf == false && t_hi < 2*n) : true);
    }

    // return the non-minimized result
    return t_hi;
}


template <typename T>
struct DefaultRedcLargeR
{
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, FullrangeTag)
  {
    bool ovf;
    T result = REDC_non_minimized(ovf, u_hi, u_lo, n, neg_inv_n);
    // REDC_non_minimized()'s Postcondition #1 guarantees the following
    T minimized_result = (ovf || result>=n) ? static_cast<T>(result-n) : result;
    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
  }

  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, HalfrangeTag)
  {
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // HalfrangeTag has the precondition requirement that n < R/2 (see
        // MontyHalfRange for more on this).
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(n < Rdiv2);
    }

    bool ovf;
    T result = REDC_non_minimized(ovf, u_hi, u_lo, n, neg_inv_n);
    // Since we have the precondition n<R/2, we know from REDC_non_minimized()'s
    // Postcondition #4 that ovf is false.
    HPBC_ASSERT2(ovf == false);
    // Since ovf == false, REDC_non_minimized()'s Postcondition #1 guarantees
    T minimized_result = (result >= n) ? static_cast<T>(result - n) : result;
    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
  }

  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, QuarterrangeTag)
  {
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // QuarterrangeTag has the precondition requirement that n < R/4 (see
        // MontyQuarterRange for more on this).
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(n < Rdiv4);
    }

    bool ovf;
    T result = REDC_non_minimized(ovf, u_hi, u_lo, n, neg_inv_n);
    // Since we have the precondition n<R/4, we know from REDC_non_minimized()'s
    // Postcondition #4 that ovf is false and result < 2*n.
    HPBC_ASSERT2(ovf == false);
    HPBC_POSTCONDITION2(result < 2*n);
    // MontyQuarterRange (and hence QuarterrangeTag) allows any montgomery
    // values that satisfy 0 <= value < 2*n, so this result doesn't need to be
    // further reduced.
    return result;
  }


  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, FullrangeTag)
  {
    // any MontyFullRange (and MontyHalfRange) montgomery value should be < n
    HPBC_PRECONDITION2(x < n);

    T u_hi = 0;  T u_lo = x;
    bool ovf;
    T thi = REDC_non_minimized(ovf, u_hi, u_lo, n, neg_inv_n);
    // Since we have u_hi == 0 and u_lo < n, we know from REDC_non_minimized()'s
    // Postcondition #2 that ovf == false, and from Postcondition #3 that
    // thi < n.
    HPBC_ASSERT2(ovf == false && thi < n);
    // Combining this with REDC_non_minimized() Postcondition #1 gives us
    T minimized_result = thi;

    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
  }

  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, HalfrangeTag)
  {
    // The implementations for Halfrange and Fullrange should be the exact same.
    return convert_out(x, n, neg_inv_n, FullrangeTag());
  }
  
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, QuarterrangeTag)
  {
    // QuarterrangeTag implies a requirement that n < R/4, but this function
    // will work for any n < R (since n is a T value, it satisfies this).  Since
    // we don't actually need n < R/4, we don't make it a precondition.

    T u_hi = 0;  T u_lo = x;
    bool ovf;
    T thi = REDC_non_minimized(ovf, u_hi, u_lo, n, neg_inv_n);
    // Since we have u_hi == 0, we know from REDC_non_minimized()'s
    // Postcondition #2 that ovf is false.
    HPBC_ASSERT2(ovf == false);
    // Given that ovf==false, REDC_non_minimized()'s Postcondition #1 guarantees
    T minimized_result = (thi >= n) ? (thi - n) : thi;

    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
  }
};

} // end namespace detail_redc_large



// ----------------
// Public Functions
// ----------------


// primary template
template <typename T>
struct RedcLargeR
{
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
  static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");

  // MTAG (montgomery tag) can be FullrangeTag, HalfrangeTag, QuarterrangeTag
  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, MTAG, PTAG)
  {
    return detail_redc_large::
                   DefaultRedcLargeR<T>::REDC(u_hi, u_lo, n, neg_inv_n, MTAG());
  }

  template <class MTAG>
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, MTAG)
  {
    return detail_redc_large::
                     DefaultRedcLargeR<T>::convert_out(x, n, neg_inv_n, MTAG());
  }
};



#if defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

// specialization for uint64_t (for x86_64)
template <>
struct RedcLargeR<std::uint64_t>
{
  using T = std::uint64_t;

  // This version should have: cycles latency 11, fused uops 12
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, FullrangeTag, InplaceLowlatencyTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...FullrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T rrdx, dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
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
    T result = rrax + rrdx;   // let compiler choose between add/lea
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                     u_hi, u_lo, n, neg_inv_n, FullrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }

  // This version should have: cycles latency 11, fused uops 11
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, FullrangeTag,OutofplaceLowlatencyTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...FullrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T reg = u_hi;
    T rrdx, dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
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
    T result = reg;
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                     u_hi, u_lo, n, neg_inv_n, FullrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }

  // This next version is expected to have:
  // cycles latency a little over 11 (perhaps ~11.2 on average), fused uops 11.
  // Note:
  // The expected latency here is tricky, and depends on mul.  The uops.info
  // latency figures for the mul instruction suggest that rax will be written by
  // a uop (as part of the mul) one cycle earlier than a separate uop which
  // writes rdx - so that rax takes 3 cycles to produce and rdx takes 4 cycles.
  // My testing supports this (I used an asm loop containing mul, with each loop
  // iteration dependent upon either rax or rdx produced by the previous mul).
  // Assuming this is indeed true, we should expect a latency of 11 for this
  // REDC.  However, I consistently measure the corresponding low-latency
  // version of this function as being roughly 1.5% faster than this (low-uops)
  // version, using a latency bound loop test on a Haswell CPU, despite the fact
  // that we would have expected both functions to have the same latency of 11
  // cycles, and despite both functions containing close to identical assembly
  // code.  Although such a small difference may often be explained by quirks
  // in the surrounding assembly, and perhaps alignment, this appears to be a
  // real performance difference - I see it with surrounding assembly that looks
  // exactly the same (though of course the contents of the inlined asm function
  // differ), and with gcc and clang, and when switching the order of the test
  // loops.  A possible explanation might be found in looking at this function's
  // asm:  It has a very tight 1 cycle window from when rax gets produced by the
  // mul, during which the addq can and must execute, in order for both the addq
  // result and rdx (from the mul) to both be ready at the same time to then be
  // used by the adcq instruction.  If on occasion addq doesn't execute in the
  // same cycle in which rax becomes available, on average a fraction of one
  // cycle of latency would be incurred.  For comparison, the corresponding low
  // latency version of this function uses xorlq and negq instead of addq, and
  // it has a roughly 6 cycle window in which it can execute those two
  // instructions (without incurring extra latency).  This is a guess, and the
  // true reason may be different - perf or vtune would help.
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, FullrangeTag, InplaceLowuopsTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...FullrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T rrdx, dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
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
    T result = rrax + rrdx;   // let compiler choose between add/lea
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                     u_hi, u_lo, n, neg_inv_n, FullrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }

  // This version is expected to have:
  // cycles latency a little over 11 (~11.2 on average), fused uops 10.
  // See the comment for the FullrangeTag/InplaceLowuopsTag version of this
  // function above, for details on why expected latency has a fractional value.
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, FullrangeTag, OutofplaceLowuopsTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...FullrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T reg = u_hi;
    T rrdx, dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
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
    T result = reg;
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                     u_hi, u_lo, n, neg_inv_n, FullrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }



  // Note: PrivateInplaceTag covers InplaceLowlatencyTag and InplaceLowuopsTag.
  // This version should have: cycles latency 10, fused uops 9
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, HalfrangeTag, PrivateInplaceTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // HalfrangeTag implicitly has the precondition requirement that
        // n < R/2 (see MontyHalfRange for the explicit precondition).
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(n < Rdiv2);
    }

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...HalfrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
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
    T result = rrax;
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                     u_hi, u_lo, n, neg_inv_n, HalfrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }

  // PrivateOutofplaceTag covers OutofplaceLowlatencyTag / OutofplaceLowuopsTag.
  // This version should have: cycles latency 10, fused uops 9
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, HalfrangeTag, PrivateOutofplaceTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // HalfrangeTag implicitly has the precondition requirement that
        // n < R/2 (see MontyHalfRange for the explicit precondition).
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(n < Rdiv2);
    }

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...HalfrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T reg = u_hi;
    T dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
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
    T result = reg;
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                     u_hi, u_lo, n, neg_inv_n, HalfrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }



  // Note: PrivateInplaceTag covers InplaceLowlatencyTag and InplaceLowuopsTag.
  // This version should have: cycles latency 8, fused uops 7
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, QuarterrangeTag, PrivateInplaceTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // QuarterrangeTag implicitly has the precondition requirement that
        // n < R/4 (see MontyQuarterRange for the explicit precondition).
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(n < Rdiv4);
    }

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...QuarterrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, sum==0 */
          "movq %[uhi], %%rax \n\t" /* rax = u_hi */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %%rdx, %%rax \n\t"    /* t_hi = addcarry(u_hi, mn_hi) */
        : "+&a"(rrax), [tmp]"=&r"(dummy)
        : [uhi]"r"(u_hi), [n]"rm"(n), [inv]"rm"(neg_inv_n)
        : "rdx", "cc");
    T result = rrax;
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                  u_hi, u_lo, n, neg_inv_n, QuarterrangeTag()));

    HPBC_POSTCONDITION2(result < 2*n);
    return result;
  }

  // PrivateOutofplaceTag covers OutofplaceLowlatencyTag / OutofplaceLowuopsTag.
  // This version should have: cycles latency 8, fused uops 6
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T neg_inv_n, QuarterrangeTag,PrivateOutofplaceTag)
  {
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_minimized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // QuarterrangeTag implicitly has the precondition requirement that
        // n < R/4 (see MontyQuarterRange for the explicit precondition).
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(n < Rdiv4);
    }

    // We'll ignore this function's calling convention since it's inlined.
    // Assume u_lo is in register rax.  Note: a 128 bit mul normally precedes
    // this function, putting u_lo in rax.
    // This asm implements DefaultRedcLargeR<uint64_t>::REDC(...QuarterrangeTag)
    // Thus, the algorithm should be correct for the same reasons given there.
    T rrax = u_lo;
    T reg = u_hi;
    T dummy;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"   /* tmp = u_lo */
        "imulq %[inv], %%rax \n\t"  /* m = u_lo * neg_inv_n */
        "mulq %[n] \n\t"            /* mn_hilo = m * n */
        /* "addq %[tmp], %%rax \n\t" */ /* rax=u_lo+mn_lo. Sets carry, sum==0 */
          "negq %[tmp] \n\t"        /* Sets carry to (u_lo != 0) */
        "adcq %%rdx, %[reg] \n\t"   /* t_hi = addcarry(u_hi, mn_hi) */
        : [reg]"+&r"(reg), "+&a"(rrax), [tmp]"=&r"(dummy)
        : [n]"rm"(n), [inv]"rm"(neg_inv_n)
        : "rdx", "cc");
    T result = reg;
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                  u_hi, u_lo, n, neg_inv_n, QuarterrangeTag()));

    HPBC_POSTCONDITION2(result < 2*n);
    return result;
  }



  // For now, I don't plan to write an x86_64 asm version for convert_out().
  // Compilers seem to generate ok code from the default C++ implementation, and
  // this is also probably unlikely to be used in performance critical loops.
  template <class MTAG>
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, MTAG)
  {
    T result = detail_redc_large::
                     DefaultRedcLargeR<T>::convert_out(x, n, neg_inv_n, MTAG());
    HPBC_POSTCONDITION2(result < n);
    return result;
  }
};

#endif   // defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) &&
         // defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)


}} // end namespace

#endif
