// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_REDC_LARGE_R_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_REDC_LARGE_R_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace montgomery_arithmetic {


namespace detail_redc_large {
// -----------------
// Private Functions
// -----------------

#if 0
// The old REDC function, using the traditional/conventional REDC algorithm
// (utilizing the negative inverse of n, mod R) from Peter Montgomery's 1985
// paper.  You can find proofs for its postconditions in older commits of this
// file in git, where this function is not disabled with #if 0.
// This function/algorithm is superceded by the better function/algo below it.
template <typename T> HURCHALLA_FORCE_INLINE
T REDC_non_minimized(bool& ovf, T u_hi, T u_lo, T n, T neg_inv_n)
{
    T m = u_lo * neg_inv_n;  // computes  m = (u * neg_inv_n) % R
    T mn_lo;
    T mn_hi = unsigned_multiply_to_hilo_product(&mn_lo, m, n);
    T t_hi = u_hi + mn_hi;
    t_hi += (u_lo != 0);
    ovf = (t_hi < u_hi);

    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (ovf || t_hi >= n) ? (t_hi - n) : t_hi;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (ovf == false && t_hi < 2*n) : true);
    }
    HPBC_POSTCONDITION2((u_hi == 0) ? ovf == false : true);
    HPBC_POSTCONDITION2((u_hi == 0 && u_lo < n) ? t_hi < n : true);
    // I strongly suspect I can prove a related postcondition, but I'm not
    // aware of any additional benefit it would provide-
    //    HPBC_POSTCONDITION2((u_hi == 0) ? t_hi <= n : true);
    return t_hi;
}
#endif


// Generic Montgomery REDC algorithm

// This function implements the REDC algorithm as described at
// hurchalla/montgomery_arithmetic/detail/README_REDC.md.
// This is an alternate version of the REDC algorithm, that differs in small but
// important ways from Peter Montgomery's original 1985 paper "Modular
// multiplication without trial division".  For our purposes, the most important
// distinction is that it is a more efficient algorithm both for latency and
// number of instructions.  See README_REDC.md for more info.
// Note that the description in README_REDC.md uses a variable name "T", despite
// the fact that in C++ "T" is conventionally reserved for use as a template
// parameter name.   This is done for consistency with most all descriptions of
// Montgomery multiplication/REDC, including Montgomery's 1985 paper, the
// Wikipedia article https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
// and many more.
//
// In the function below, we will use the variable name "u" (and u_hi and u_lo)
// in place of the algorithm description's variable name "T", since as already
// stated, in C++ "T" is de-facto a reserved name for a template parameter.  We
// will use "T" in its normal C++ sense, as a template parameter name - this has
// no relation to "T" in the algorithm description.
// We also use "n" instead of "N", and "inv_n" instead of "N^(-1)" (N with the
// superscript -1).  The constant "R" remains the same, and represents the value
// R = 2^(ma::ma_numeric_limits<T>::digits).  As an example, if T is uint64_t,
// then R = 2^64.
// This function's name is "REDC_non_finalized" to reflect the fact that the
// value it returns is not finalized to the least residual class mod n (i.e. so
// that 0 <= return_value < n).  See the starred* comment under Postcondition #1
// for more info.
template <typename T> HURCHALLA_FORCE_INLINE
T REDC_non_finalized(bool& ovf, T u_hi, T u_lo, T n, T inv_n)
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

    // assert(n * inv_n ≡ 1 (mod R))
    HPBC_PRECONDITION2(static_cast<T>(n * inv_n) == 1);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T m = u_lo * inv_n;  // computes  m = (u * inv_n) % R

    T mn_lo;
    T mn_hi = unsigned_multiply_to_hilo_product(&mn_lo, m, n);

    // mn = m*n.  Since m = (u_lo*inv_n)%R, we know m < R, and thus  mn < R*n.
    // Therefore mn == mn_hi*R + mn_lo < R*n, and mn_hi*R < R*n - mn_lo <= R*n,
    // and thus  mn_hi < n.
        // *** Assertion #1 ***
    HPBC_ASSERT2(mn_hi < n);

    // Compute (u - mn)/R :
    T t_hi = u_hi - mn_hi;   // computes  t_hi = (u_hi - mn_hi) % R

    ovf = (u_hi < mn_hi);    // lets us know if the subtraction overflowed

    // We do not need to explicitly perform the low part subtraction
    // (u_lo - mn_lo), because the REDC algorithm guarantees
    // (u_lo - mn_lo) % R == 0.  Since both u_lo < R and mn_lo < R, this means
    // that u_lo == mn_lo, and thus (u_lo - mn_lo) will never generate a borrow/
    // carry.  We will simply ignore this low part subtraction.
        // *** Assertion #2 ***
    HPBC_ASSERT2(u_lo == mn_lo);

    // Since u_hi and u_lo are type T (which is unsigned) variables, both
    // u_hi >= 0 and u_lo >= 0, and thus  u = u_hi*R + u_lo >= 0.  Along with
    // the precondition u < n*R, we therefore know that  0 <= u < n*R.
    // Since m is type T we know m >= 0, and since n is type T we know n >= 0.
    // Thus mn >= 0.  Subtracting mn from all parts of 0 <= u < n*R, we have
    // -mn <= u - mn < n*R - mn.  Since we just showed mn >= 0, we know
    // -mn <= 0 and thus  n*R - mn <= n*R.  Assertion #1 states that mn < n*R,
    // and thus we know  -n*R < -mn.  Therefore
    // -n*R < -mn <= u - mn < n*R - mn <= n*R.  Simplified, -n*R < u - mn < n*R.
    // Since (u - mn) is divisible by R (see README_REDC.md for proof), and
    // R > 0, we have  -n < (u - mn)/R < n.
    // Since u = u_hi*R + u_lo, and mn = mn_hi*R + mn_lo,  we know
    // u - mn == u_hi*R + u_lo - mn_hi*R - mn_lo.  Assertion #2 states that
    // u_lo == mn_lo, and so we have
    // u - mn == u_hi*R - mn_hi*R == (u_hi - mn_hi)*R,  and thus
    // (u - mn)/R == u_hi - mn_hi.  Therefore,  -n < u_hi - mn_hi < n.
    //
    // All this translates into
    // Postcondition #1
    // ----------------
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T finalized_result = (ovf) ? static_cast<T>(t_hi + n) : t_hi;
        HPBC_POSTCONDITION2(finalized_result < n);
    }
    // * Aside from this postcondition, we do not actually compute the finalized
    // least residual mod n result, because some Montgomery Forms are
    // constrained in ways that allow a simpler and more efficient computation
    // of the finalized result.  For example, in some forms the input u_hi (and
    // the return value) is allowed to occupy the range 0 <= u_hi < 2*n, which
    // lets us change the conditional add of n at the end of the REDC into an
    // unconditional add of n.

    // Postcondition #2:  If  n < R/2,  then  0 < t_hi + n < 2*n
    // ---------------------------------------------------------
    // We already showed  -n < u_hi - mn_hi < n.  Adding n to all parts we get
    // 0 < u_hi - mn_hi + n < 2*n.  Although this is true regardless of the size
    // of n, we can only test this postcondition when n < R/2 (any larger value
    // of n would overflow on 2*n).
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (0 < static_cast<T>(t_hi + n)) &&
                                       (static_cast<T>(t_hi + n) < 2*n) : true);
    }

    // return the non-finalized result
    return t_hi;
}



template <typename T>
struct DefaultRedcLargeR
{
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T inv_n, FullrangeTag)
  {
    bool ovf;
    T result = REDC_non_finalized(ovf, u_hi, u_lo, n, inv_n);
    // REDC_non_finalized()'s Postcondition #1 guarantees the following
    T finalized_result = (ovf) ? static_cast<T>(result + n) : result;
    HPBC_POSTCONDITION2(finalized_result < n);
    return finalized_result;
  }
  // Since HalfrangeTag inherits from FullrangeTag, it will argument match to
  // the FullrangeTag function above.  I don't see any way to improve upon
  // the FullrangeTag function with a dedicated version for HalfrangeTag.


  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T inv_n, QuarterrangeTag)
  {
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // QuarterrangeTag has the precondition requirement that n < R/4 (see
        // MontyQuarterRange for more on this).
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(n < Rdiv4);
    }

    bool ovf;
    T result= static_cast<T>(n + REDC_non_finalized(ovf, u_hi, u_lo, n, inv_n));
    // Since we have the precondition n<R/4, we know from REDC_non_finalized()'s
    // Postcondition #2 that  0 < result < 2*n.
    HPBC_POSTCONDITION2(0 < result && result < 2*n);
    // MontyQuarterRange (and hence QuarterrangeTag) allows any montgomery
    // values that satisfy 0 <= value < 2*n, so this result doesn't need to be
    // further reduced.
    return result;
  }
  // Since SixthrangeTag inherits from QuarterrangeTag, it will argument match
  // to the QuarterrangeTag function above.


  // Note that the old REDC_non_finalized function (disabled at top) can provide
  // superior performance for use with convert_out, for FullrangeTag and
  // HalfrangeTag (I don't believe it would benefit QuarterrangeTag or
  // SixthrangeTag).  If we decided at some point to change convert_out to use
  // the old REDC_non_finalized, we would need to call the old REDC with -inv_n,
  // since it requires the negative inverse for an argument.  The negation
  // requires an extra instruction that would take a cycle most likely - we
  // would probably only care about performance in a loop context and the
  // compiler hopefully would hoist the negation out of the loop in that case,
  // but in the worst case scenario it's possible the extra instruction for
  // negation could cancel out the benefit from using the old REDC for
  // convert_out().
  // Since convert_out() is typically not performance critical, for the moment I
  // I am electing to keep things simple by not re-introducing/re-enabling the
  // old REDC_non_finalized function.  Using the old REDC_non_finalized would
  // likely save 1 or 2 cycles, but convert_out() doesn't seem important enough
  // to justify it.
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T inv_n)
  {
    T u_hi = 0;  T u_lo = x;
    bool ovf;
    T result = REDC_non_finalized(ovf, u_hi, u_lo, n, inv_n);
    // REDC_non_finalized()'s Postcondition #1 guarantees the following
    T finalized_result = (ovf) ? static_cast<T>(result + n) : result;
    HPBC_POSTCONDITION2(finalized_result < n);
    return finalized_result;
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

  // For MTAGs see monty_tag_structs.h; for PTAGs see optimization_tag_structs.h
  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T inv_n, MTAG, PTAG)
  {
    return detail_redc_large::
                       DefaultRedcLargeR<T>::REDC(u_hi, u_lo, n, inv_n, MTAG());
  }

  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T inv_n)
  {
    return detail_redc_large::DefaultRedcLargeR<T>::convert_out(x, n, inv_n);
  }
};



#if defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

// specialization for uint64_t (for x86_64)
template <>
struct RedcLargeR<std::uint64_t>
{
  using T = std::uint64_t;

  // This version should have: cycles latency 9, fused uops 7
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T inv_n, FullrangeTag, LowlatencyTag)
  {
    // This implementation is based closely on DefaultRedcLargeR<uint64_t>::REDC
    // for FullrangeTag.  Thus the algorithm should be correct for the same
    // reasons given there.
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_finalized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);
    HPBC_PRECONDITION2(static_cast<T>(n * inv_n) == 1);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T m = u_lo * inv_n;
    T mn_lo;
    T mn_hi = unsigned_multiply_to_hilo_product(&mn_lo, m, n);
    T reg = u_hi + n;
    T uhi = u_hi;
    __asm__ (
        "subq %[mnhi], %[reg] \n\t"     /* reg = u_hi + n - mn_hi */
        "subq %[mnhi], %[uhi] \n\t"     /* t_hi = u_hi - mn_hi */
        "cmovaeq %[uhi], %[reg] \n\t"   /* reg = (u_hi >= mn_hi) ? t_hi : reg */
        : [reg]"+&r"(reg), [uhi]"+&r"(uhi)
        : [mnhi]"r"(mn_hi)
        : "cc");
    T result = reg;
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                         u_hi, u_lo, n, inv_n, FullrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }

  // This version should have: cycles latency 10, fused uops 6
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T inv_n, FullrangeTag, LowuopsTag)
  {
    // This implementation is based closely on DefaultRedcLargeR<uint64_t>::REDC
    // for FullrangeTag.  Thus the algorithm should be correct for the same
    // reasons given there.
    // We require u = (u_hi*R + u_lo) < n*R.  As shown in precondition #1 in
    // REDC_non_finalized(), u_hi < n guarantees this.
    HPBC_PRECONDITION2(u_hi < n);
    HPBC_PRECONDITION2(static_cast<T>(n * inv_n) == 1);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T m = u_lo * inv_n;
    T mn_lo;
    T mn_hi = unsigned_multiply_to_hilo_product(&mn_lo, m, n);
    namespace ma = hurchalla::modular_arithmetic;
    T result = ma::modular_subtraction_prereduced_inputs(u_hi, mn_hi, n);
    HPBC_ASSERT2(result == detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                         u_hi, u_lo, n, inv_n, FullrangeTag()));

    HPBC_POSTCONDITION2(result < n);
    return result;
  }

  // We don't need a dedicated REDC function for HalfrangeTag since
  // FullrangeTag is already optimal for it, and HalfrangeTag will match to it.


  // This version should have: cycles latency 8, fused uops 5.
  // DefaultRedcLargeR's REDC for QuarterrangeTag should already be optimal.
  static HURCHALLA_FORCE_INLINE
  T REDC(T u_hi, T u_lo, T n, T inv_n, QuarterrangeTag, PrivateAnyTag)
  {
    return detail_redc_large::DefaultRedcLargeR<T>::REDC(
                                       u_hi, u_lo, n, inv_n, QuarterrangeTag());
  }

  // We don't need a dedicated REDC function for SixthrangeTag since
  // QuarterrangeTag is already optimal for it, and SixthrangeTag will match it.


  // For now, I have no plan to write an x86_64 asm version for convert_out().
  // Compilers seem to generate ok code from the default C++ implementation, and
  // this function is probably unlikely to be used in performance critical loops
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T inv_n)
  {
    T result= detail_redc_large::DefaultRedcLargeR<T>::convert_out(x, n, inv_n);
    HPBC_POSTCONDITION2(result < n);
    return result;
  }
};

#endif   // defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) &&
         // defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
