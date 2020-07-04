// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/montgomery_arithmetic/detail/sized_uint.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace montgomery_arithmetic {


// This file provides the following two helper functions to callers:
//
// T montmul_non_minimized(bool& ovf, T x, T y, T n, T neg_inv_n) -
//      Multiplies two mongomery values x and y, and returns their non-minimized
//      (mod n) montgomery product.
//
// T montout_non_minimized(T x, T n, T neg_inv_n) -
//      Converts the montgomery value 'x' to a non-minimized (mod n) standard
//      integer.


namespace detail_monty_common {
// -----------------
// Private Functions
// -----------------


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
template <typename T>
HURCHALLA_FORCE_INLINE
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

    // assert(n * neg_inv_n ≡ -1 (mod R))
    HPBC_PRECONDITION2(n * neg_inv_n == static_cast<T>(0) - static_cast<T>(1));
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    // We require the precondition  u < n*R.  Or  u == u_hi*R + u_lo < n*R.
    // If u_hi < n:  then u_hi+1 <= n, and u_hi*R + R <= n*R.  Since u_lo < R,
    //   u == u_hi*R + u_lo < u_hi*R + R <= n*R.  We would have  u < n*R, and
    //   so u_hi < n is always sufficient to satisfy the precondition.
    // If u_hi >= n:  then u_hi*R >= n*R, and u == u_hi*R + u_lo >= n*R, which
    //   fails the precondition.
    // Thus u_hi < n is sufficient and necessary to satisfy the precondition.
    HPBC_PRECONDITION2(u_hi < n);

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
    // We allow the Montgomery Form calling this function to compute its final
    // minimized result in a potentially more efficient manner.


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
    // [ These findings were inspired by ideas in section 5 of "Montgomery's
    // Multiplication Technique: How to Make It Smaller and Faster", at
    // https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps ]
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (ovf == false && t_hi < 2*n) : true);
    }

    // return the non-minimized result
    return t_hi;
}


// This function is intended for use when T is smaller than the ALU's native bit
//   width.  It multiplies two mongomery values x and y.  It is based on the
//   REDC algorithm (but it uses a multiplication product at the function start
//   as the input to REDC, and it does not minimize the output result).  See
//   https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
// Precondition: with theoretical unlimited precision standard multiplication,
//   we require  x*y < n*R.  The constant R represents the value
//   R = 2^(ma::ma_numeric_limits<T>::digits).
// Returns: the non-minimized montgomery product.
// Implmentation: this function is a safer and generic version of the following-
// uint32_t montmul_non_minimized(bool& ovf, uint32_t x, uint32_t y, uint32_t n,
//                                                           uint32_t neg_inv_n)
// {
//    using T = uint32_t;
//    using T2 = uint64_t;
//    T2 u = static_cast<T2>(x) * static_cast<T2>(y);
//    T m = static_cast<T>(u) * neg_inv_n;
//    T2 mn = static_cast<T2>(m) * static_cast<T2>(n);
//    T2 t = u + mn;
//    ovf = (t < u);
//    T t_hi = static_cast<T>(t >> ma::ma_numeric_limits<T>::digits);
//    return t_hi;
// }
template <typename T, typename T2>
HURCHALLA_FORCE_INLINE
T typecast_montmul_non_minimized(bool& ovf, T x, T y, T n, T neg_inv_n)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T>::is_modulo, "");

    static_assert(ma::ma_numeric_limits<T2>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T2>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T2>::is_modulo, "");

    static_assert(ma::ma_numeric_limits<T2>::digits ==
                  2 * ma::ma_numeric_limits<T>::digits, "");

    // For casts, we want to use types that are protected from surprises and
    // undefined behavior due to the unsigned integral promotion rules in C++.
    // https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/
    using V = typename safely_promote_unsigned<T>::type;
    using V2 = typename safely_promote_unsigned<T2>::type;
    static_assert(ma::ma_numeric_limits<V>::is_modulo, "");
    static_assert(ma::ma_numeric_limits<V2>::is_modulo, "");

    static constexpr int bit_width_T = ma::ma_numeric_limits<T>::digits;
    static constexpr V2 R = static_cast<V2>(1) << bit_width_T;

    // assert(n * neg_inv_n ≡ -1 (mod R))
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<V>(n) * static_cast<V>(neg_inv_n)) ==
                static_cast<T>(static_cast<V>(0) - static_cast<V>(1))
                );
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    V2 u = static_cast<V2>(x) * static_cast<V2>(y);
    HPBC_PRECONDITION2(u < static_cast<V2>(n) * R);

    T m = static_cast<T>(static_cast<V>(u) * static_cast<V>(neg_inv_n));
    V2 mn = static_cast<V2>(m) * static_cast<V2>(n);

    T2 t = static_cast<T2>(u + mn);
    ovf = (static_cast<V2>(t) < u);

    T t_hi = static_cast<T>(static_cast<V2>(t) >> bit_width_T);

    // For the same reason as REDC_non_minimized(), we do not compute the final
    // minimized result, aside from this postcondition.
    // Postcondition #5
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (ovf || t_hi >= n) ? static_cast<T>(t_hi-n) : t_hi;
        HPBC_POSTCONDITION2(minimized_result < n);
    }

    // Postcondition #6 - If y == 1, then ovf == false.
    // ------------------------------------------------
    // We know mn <= (R-1)*(R-1).  For this case, y==1, so  u = x*y == x <= R-1.
    // And so u + mn <= (R-1) + ((R-1)*(R-1)) == (R-1)*R < R*R.  Since we also
    // know  u + mn >= 0, we know  0 <= u + mn < R*R, and thus
    // (u + mn) == (u + mn) % (R*R).  Therefore since t = (u + mn) % (R*R),
    // t == u + mn, and thus t >= u.  This means ovf = (t < u) must be false.
    HPBC_POSTCONDITION2((y==1) ? ovf==false : true);

    // Postcondition #7 - If y == 1 and x < n, then t_hi < n.
    // ------------------------------------------------------
    // For this case, y==1 and x<n, and so u = x*y == x < n.  Thus u <= n-1.
    // And since m is type T, m <= R-1.  Therefore
    // u + m*n <= (n-1) + (R-1)*n == R*n - 1,  and thus  u + m*n < R*n < R*R.
    // Using similar reasoning to Postcondition #6, we therefore know
    // t = (u + mn) % (R*R) == (u + mn),  and since  u + m*n < R*n,  t < R*n.
    // Since t_hi = t/R (which divides evenly), we know
    // t_hi = t/R < (R*n)/R == n, and thus,  t_hi < n.
    // [ This finding was inspired by ideas in section 5 of "Montgomery's
    // Multiplication Technique: How to Make It Smaller and Faster", at
    // https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps ]
    HPBC_POSTCONDITION2((y==1 && x<n) ? t_hi<n : true);

    // Postcondition #8 - If n < R/2, then ovf == false  and  t_hi < 2*n.
    // ------------------------------------------------------------------
    // Since m is type T we know m < R, and since mn = m*n,  mn < R*n.  Adding
    // this with the function precondition u < n*R, we get  u + mn < 2*R*n.
    // Since this case specifies n < R/2, we have  u + mn < 2*R*R/2 == R*R.
    // Using similar reasoning to Postcondition #6, we therefore know
    // t = (u + mn) % (R*R) == (u + mn).  Thus t >= u.  Therefore  ovf = (t < u)
    // must be false.
    // We've already shown u + mn < 2*R*n, and since t == (u + mn), we therefore
    // know  t < 2*R*n.  Since t_hi = t/R (which divides evenly), we know
    // t_hi = t/R < (2*R*n)/R == 2*n, and thus,  t_hi < 2*n.
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ?
                       (ovf == false && t_hi < static_cast<T>(2)*n) : true);
    }

    // return the non-minimized result
    return t_hi;
}




// Multiplies two mongomery values x and y.  Requires x*y < n*R (assuming
// theoretical infinite precision standard multiplication).
// For postconditions, see inside this function.
// Returns the non-minimized (mod n) montgomery product.
template <typename T>
HURCHALLA_FORCE_INLINE
T impl_montmul_non_minimized(bool& ovf, T x, T y, T n, T neg_inv_n)
{
    namespace ma = hurchalla::modular_arithmetic;
    T u_lo;
    T u_hi = unsigned_multiply_to_hilo_product(&u_lo, x, y);
    // Having u_hi < n is sufficient and necessary to satisfy the requirement of
    // x*y == u < n*R.  See REDC_non_minimized() for proof.
    HPBC_PRECONDITION2(u_hi < n);

    T result = REDC_non_minimized(ovf, u_hi, u_lo, n, neg_inv_n);

    // REDC_non_minimized Postcondition #4 guarantees that if n < R/2, then
    // ovf == false  and  result < 2*n.
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (ovf==false && result<2*n) : true);
    }
    // REDC_non_minimized Postcondition #1 guarantees the following
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (ovf || result >= n) ? (result - n) : result;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    return result;
}

// Converts the montgomery value 'x' to a non-minimized (mod n) standard integer
// For postconditions, see inside this function.
template <typename T>
HURCHALLA_FORCE_INLINE T impl_montout_non_minimized(T x, T n, T neg_inv_n)
{
    T u_hi = static_cast<T>(0);
    T u_lo = x;
    bool ovf;
    // Having u_hi == 0 satisfies REDC_non_minimized's precondition that u < n*R

    T result = REDC_non_minimized(ovf, u_hi, u_lo, n, neg_inv_n);
    // REDC_non_minimized Postcondition #2 guarantees ovf==false, since u_hi==0.
    HPBC_ASSERT2(ovf == false);

    // REDC_non_minimized Postcondition #3 guarantees the following, for u_hi==0
    HPBC_POSTCONDITION2((x < n) ? result < n : true);
    // REDC_non_minimized Postcondition #1 guarantees the following, since ovf
    // is false:
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (result >= n) ? (result - n) : result;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    return result;
}

// -----------------------------------------------------------------------------
// -- Non-template overloads (they have first priority for argument matching) --
// -----------------------------------------------------------------------------

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
#endif

#if HURCHALLA_TARGET_BIT_WIDTH >= 16
HURCHALLA_FORCE_INLINE uint8_t impl_montmul_non_minimized(bool& ovf, uint8_t x,
                                        uint8_t y, uint8_t n, uint8_t neg_inv_n)
{
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  x*y < n*R.
    namespace ma = hurchalla::modular_arithmetic;
    using T = uint8_t;
    using T2 = uint16_t;
    T result = typecast_montmul_non_minimized<T, T2>(ovf, x, y, n, neg_inv_n);

    // typecast_montmul_non_minimized Postcondition #8 guarantees that
    // if n < R/2, then  ovf == false  and  result < 2*n.
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (ovf==false && result<2*n) : true);
    }
    // typecast_montmul_non_minimized Postcondition #5 guarantees the following
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (ovf || result >= n) ? static_cast<T>(result - n) :
                                                    result;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    // Note that these postconditions are the same as the template function
    // version of impl_montmul_non_minimized.
    return result;
}
HURCHALLA_FORCE_INLINE
uint8_t impl_montout_non_minimized(uint8_t x, uint8_t n, uint8_t neg_inv_n)
{
    using T = uint8_t;
    using T2 = uint16_t;
    bool ovf;
    T y = static_cast<T>(1);
    // Having y == 1 satisfies typecast_montmul_non_minimized's precondition
    // that x*y < n*R, since x < R < n*R.

    T result = typecast_montmul_non_minimized<T, T2>(ovf, x, y, n, neg_inv_n);
    // typecast_montmul_non_minimized Postcondition #6 guarantees that since
    // y==1, ovf==false.
    HPBC_ASSERT2(ovf == false);

    // typecast_montmul_non_minimized Postcondition #7 guarantees the following,
    // since y==1:
    HPBC_POSTCONDITION2((x < n) ? result < n : true);
    // typecast_montmul_non_minimized Postcondition #5 guarantees the following,
    // since we know ovf == false:
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (result >= n) ? static_cast<T>(result-n) : result;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    // Note that these postconditions are the same as the template function
    // version of impl_montout_non_minimized.
    return result;
}
#endif

#if HURCHALLA_TARGET_BIT_WIDTH >= 32
HURCHALLA_FORCE_INLINE uint16_t impl_montmul_non_minimized(bool& ovf,
                         uint16_t x, uint16_t y, uint16_t n, uint16_t neg_inv_n)
{
    using T = uint16_t;
    using T2 = uint32_t;
    T result = typecast_montmul_non_minimized<T, T2>(ovf, x, y, n, neg_inv_n);
    // Postconditions are the same as for the uint8_t overload of this function.
    return result;
}
HURCHALLA_FORCE_INLINE
uint16_t impl_montout_non_minimized(uint16_t x, uint16_t n, uint16_t neg_inv_n)
{
    using T = uint16_t;
    using T2 = uint32_t;
    bool ovf;
    T y = static_cast<T>(1);
    T result = typecast_montmul_non_minimized<T, T2>(ovf, x, y, n, neg_inv_n);
    HPBC_ASSERT2(ovf == false);
    // Postconditions are the same as for the uint8_t overload of this function.
    return result;
}
#endif

#if HURCHALLA_TARGET_BIT_WIDTH >= 64
HURCHALLA_FORCE_INLINE uint32_t impl_montmul_non_minimized(bool& ovf,
                         uint32_t x, uint32_t y, uint32_t n, uint32_t neg_inv_n)
{
    using T = uint32_t;
    using T2 = uint64_t;
    T result = typecast_montmul_non_minimized<T, T2>(ovf, x, y, n, neg_inv_n);
    // Postconditions are the same as for the uint8_t overload of this function.
    return result;
}
HURCHALLA_FORCE_INLINE
uint32_t impl_montout_non_minimized(uint32_t x, uint32_t n, uint32_t neg_inv_n)
{
    using T = uint32_t;
    using T2 = uint64_t;
    bool ovf;
    T y = static_cast<T>(1);
    T result = typecast_montmul_non_minimized<T, T2>(ovf, x, y, n, neg_inv_n);
    HPBC_ASSERT2(ovf == false);
    // Postconditions are the same as for the uint8_t overload of this function.
    return result;
}
#endif


} // end namespace detail_monty_common



// ----------------
// Public Functions
// ----------------


// Multiplies two mongomery values x and y, and returns their non-minimized
// (mod n) montgomery product.
template <typename T>
HURCHALLA_FORCE_INLINE
T montmul_non_minimized(bool& ovf, T x, T y, T n, T neg_inv_n)
{
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  x*y < n*R.  (The constant R
    // represents the value R = 2^(ma::ma_numeric_limits<T>::digits))

    namespace ma = hurchalla::modular_arithmetic;
    namespace dmc = detail_monty_common;
    T result = dmc::impl_montmul_non_minimized(ovf, x, y, n, neg_inv_n);

    // (All versions of impl_montmul_non_minimized guarantee the following)
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (ovf==false && result<2*n) : true);
    }
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (ovf || result >= n) ? static_cast<T>(result - n)
                                                  : result;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    return result;
}

// Converts the montgomery value 'x' to a non-minimized (mod n) standard integer
template <typename T>
HURCHALLA_FORCE_INLINE T montout_non_minimized(T x, T n, T neg_inv_n)
{
    namespace dmc = detail_monty_common;
    T result = dmc::impl_montout_non_minimized(x, n, neg_inv_n);

    // (All versions of impl_montout_non_minimized guarantee the following)
    HPBC_POSTCONDITION2((x < n) ? result < n : true);
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (result >= n) ? static_cast<T>(result-n) : result;
        HPBC_POSTCONDITION2(minimized_result < n);
    }
    return result;
}


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
