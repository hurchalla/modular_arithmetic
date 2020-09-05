// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_REDC_SMALL_R_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_REDC_SMALL_R_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/montgomery_arithmetic/detail/sized_uint.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


namespace detail_redc_small {
// -----------------
// Private Functions
// -----------------


// This function is intended for use when T is smaller than the ALU's native bit
// width.  It corresponds precisely with RedcLargeR.h's function
// detail_redc_large::REDC_non_minimized() [which is intended for T greater than
// or equal to the ALU native bit width].  This function is based on the REDC
// algorithm, but it does not minimize the output result.  See
//   https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
// Precondition: with theoretical unlimited precision standard multiplication,
//   we require  u < n*R.  The constant R represents the value
//   R = 2^(ma::ma_numeric_limits<T>::digits).
// Returns: the non-minimized montgomery product.
// Implementation: this function is a safer, generic version of the following-
// uint32_t redc_non_minimized(bool& ovf, uint64_t u, uint32_t n,
//                                                           uint32_t neg_inv_n)
// {
//    using T = uint32_t;
//    using T2 = uint64_t;
//    T m = static_cast<T>(u) * neg_inv_n;
//    T2 mn = static_cast<T2>(m) * static_cast<T2>(n);
//    T2 t = u + mn;
//    ovf = (t < u);
//    T t_hi = static_cast<T>(t >> ma::ma_numeric_limits<T>::digits);
//    return t_hi;
// }
//
// Implementation note:
// ulo_input is the low half of the bits of u_input.  It may seem strange that
// the ulo_input parameter exists, since it's a duplicate of information that is
// in parameter u_input, but it allows instruction level parallelism -
// Prior operations on the high bits of u_input can execute at the same time as
// the first two multiplications in this function, since the mults depend only
// on the separate parameter ulo_input.  It would probably be clearer if this
// function took a uhi_input parameter instead of the u_input parameter, but
// this function would just need to reassemble the entire u_input anyway, since
// it never uses the high bits by themselves.
// This instruction level parallelism scenario occurs with fmsub() (which calls
// this function), and with fmadd().
template <typename T, typename T2> HURCHALLA_FORCE_INLINE
T REDC_non_minimized2(bool& ovf, T2 u_input, T ulo_input, T n, T neg_inv_n)
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

    V2 u = static_cast<V2>(u_input);
    HPBC_PRECONDITION2(u < static_cast<V2>(n) * R);

    T m = static_cast<T>(static_cast<V>(ulo_input) * static_cast<V>(neg_inv_n));
    V2 mn = static_cast<V2>(m) * static_cast<V2>(n);

    T2 t = static_cast<T2>(u + mn);
    ovf = (static_cast<V2>(t) < u);

    T t_hi = static_cast<T>(static_cast<V2>(t) >> bit_width_T);

    // For the same reasons given in Postcondition #1 in RedcLargeR.h's function
    // detail_redc_large::REDC_non_minimized(), we do not compute the final
    // minimized result, aside from here in this postcondition.
    // Postcondition #1
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (ovf || t_hi >= n) ? static_cast<T>(t_hi-n) : t_hi;
        HPBC_POSTCONDITION2(minimized_result < n);
    }

    // Postcondition #2 - If u < R, then ovf == false.
    // ---------------------------------------------------
    // We know mn <= (R-1)*(R-1).  For this case, u < R, so  u <= R-1.
    // And so u + mn <= (R-1) + ((R-1)*(R-1)) == (R-1)*R < R*R.  Since we also
    // know  u + mn >= 0, we know  0 <= u + mn < R*R, and thus
    // (u + mn) == (u + mn) % (R*R).  Therefore since t = (u + mn) % (R*R),
    // t == u + mn, and thus t >= u.  This means ovf = (t < u) must be false.
    HPBC_POSTCONDITION2((u < R) ? ovf == false : true);

    // Postcondition #3 - If u < n, then t_hi < n.
    // ------------------------------------------------------
    // For this case, u < n, so  u <= n-1.  Since m is type T,  m <= R-1.  Thus
    // u + m*n <= (n-1) + (R-1)*n == R*n - 1,  and thus  u + m*n < R*n < R*R.
    // Using similar reasoning to Postcondition #2, we therefore know
    // t = (u + mn) % (R*R) == (u + mn),  and since  u + m*n < R*n,  t < R*n.
    // Since t_hi = t/R (which divides evenly), we know
    // t_hi = t/R < (R*n)/R == n, and thus,  t_hi < n.
    HPBC_POSTCONDITION2((u < n) ? t_hi < n : true);

    // Postcondition #4 - If n < R/2, then ovf == false  and  t_hi < 2*n.
    // ------------------------------------------------------------------
    // Since m is type T we know m < R, and since mn = m*n,  mn < R*n.  Adding
    // this with the function precondition u < n*R, we get  u + mn < 2*R*n.
    // Since this case specifies n < R/2, we have  u + mn < 2*R*R/2 == R*R.
    // Using similar reasoning to Postcondition #2, we therefore know
    // t = (u + mn) % (R*R) == (u + mn).  Thus t >= u.  Therefore  ovf = (t < u)
    // must be false.
    // We've already shown u + mn < 2*R*n, and since t == (u + mn), we therefore
    // know  t < 2*R*n.  Since t_hi = t/R (which divides evenly), we know
    // t_hi = t/R < (2*R*n)/R == 2*n, and thus,  t_hi < 2*n.
    // [ This finding was inspired by ideas in section 5 of "Montgomery's
    // Multiplication Technique: How to Make It Smaller and Faster", at
    // https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps ]
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ma::ma_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ?
                       (ovf == false && t_hi < static_cast<T>(2)*n) : true);
    }

    // return the non-minimized result
    return t_hi;
}

} // end namespace detail_redc_small



// ----------------
// Public Functions
// ----------------


// Primary template.  This could be specialized for platform dependent asm code,
// though it's probably not worth the trouble.
template <typename T>
struct RedcSmallR
{
  using T2 = typename sized_uint<
                    2 * modular_arithmetic::ma_numeric_limits<T>::digits>::type;
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
  static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
  static_assert(modular_arithmetic::ma_numeric_limits<T2>::is_integer, "");
  static_assert(!(modular_arithmetic::ma_numeric_limits<T2>::is_signed), "");
  static_assert(modular_arithmetic::ma_numeric_limits<T2>::is_modulo, "");

  // Regarding function parameters u and u_lo, see the implementation note for
  // REDC_non_minimized2().

  static HURCHALLA_FORCE_INLINE
  T REDC(T2 u, T u_lo, T n, T neg_inv_n, FullrangeTag, PrivateAnyTag)
  {
    bool ovf;
    using namespace detail_redc_small;
    T thi = REDC_non_minimized2<T,T2>(ovf, u, u_lo, n, neg_inv_n);
    // We know the following from REDC_non_minimized2()'s Postcondition #1
    T minimized_result = (ovf || thi >= n) ? static_cast<T>(thi - n) : thi;
    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
  }

  static HURCHALLA_FORCE_INLINE
  T REDC(T2 u, T u_lo, T n, T neg_inv_n, HalfrangeTag, PrivateAnyTag)
  {
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // HalfrangeTag has the precondition requirement that n < R/2 (see
        // MontyHalfRange for more on this).
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(n < Rdiv2);
    }
    bool ovf;
    using namespace detail_redc_small;
    T thi = REDC_non_minimized2<T,T2>(ovf, u, u_lo, n, neg_inv_n);
    // Since we have the precondition n < R/2, we know from
    // REDC_non_minimized2()'s Postcondition #4 that ovf is false.
    HPBC_ASSERT2(ovf == false);
    // Since ovf == false, REDC_non_minimized2()'s Postcondition #1 guarantees
    T minimized_result = (thi >= n) ? static_cast<T>(thi - n) : thi;
    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
  }

  static HURCHALLA_FORCE_INLINE
  T REDC(T2 u, T u_lo, T n, T neg_inv_n, QuarterrangeTag, PrivateAnyTag)
  {
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // QuarterrangeTag has the precondition requirement that n < R/4 (see
        // MontyQuarterRange for more on this).
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(n < Rdiv4);
    }
    bool ovf;
    using namespace detail_redc_small;
    T non_min_result = REDC_non_minimized2<T,T2>(ovf, u, u_lo, n, neg_inv_n);
    // Since we have the precondition n < R/4, we know from
    // REDC_non_minimized2's Postcondition #4 that ovf is false and result < 2*n
    HPBC_ASSERT2(ovf == false);
    HPBC_POSTCONDITION2(non_min_result < 2*n);
    // MontyQuarterRange (and hence QuarterrangeTag) allows any montgomery
    // values that satisfy 0 <= value < 2*n, so this result doesn't need to be
    // further reduced.
    return non_min_result;
  }


  // Converts the montgomery value 'x' to a minimized (mod n) standard integer
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, FullrangeTag)
  {
    // MontyFullRange (and MontyHalfRange) values are always < n
    HPBC_PRECONDITION2(x < n);
    T2 u = static_cast<T2>(x);
    // u = x satisfies REDC_non_minimized2's precondition requiring u < n*R,
    // since u = x < n < R < n*R.
    bool ovf;
    using namespace detail_redc_small;
    T thi = REDC_non_minimized2<T,T2>(ovf, u, static_cast<T>(u), n, neg_inv_n);
    // Since we have u < R and u < n, we know from REDC_non_minimized2's
    // Postcondition #2 that ovf == false, and from Postcondition #3 that
    // thi < n.
    HPBC_ASSERT2(ovf == false && thi < n);
    // Combining this with REDC_non_minimized2's Postcondition #1 gives us
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
    T2 u = static_cast<T2>(x);
    // u = x satisfies REDC_non_minimized2's precondition requiring u < n*R,
    // since u = x < R < n*R.
    bool ovf;
    using namespace detail_redc_small;
    T thi = REDC_non_minimized2<T,T2>(ovf, u, static_cast<T>(u), n, neg_inv_n);
    // Since we have u < R, we know from REDC_non_minimized2's
    // Postcondition #2 that ovf is false.
    HPBC_ASSERT2(ovf == false);
    // Given that ovf==false, REDC_non_minimized2's Postcondition #1 guarantees
    T minimized_result = (thi >= n) ? static_cast<T>(thi - n) : thi;
    HPBC_POSTCONDITION2(minimized_result < n);
    return minimized_result;
  }
};


}} // end namespace


#endif
