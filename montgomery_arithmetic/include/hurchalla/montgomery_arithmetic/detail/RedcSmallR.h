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

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#  pragma warning(disable : 4309)
#endif

namespace hurchalla { namespace montgomery_arithmetic {

namespace detail_redc_small {
// -----------------
// Private Functions
// -----------------


// This function is intended for use when T is smaller than the ALU's native bit
// width.  It corresponds precisely to RedcLargeR.h's function
// detail_redc_large::REDC_non_finalized() [which is intended for T greater than
// or equal to the ALU native bit width].  This function is based on the REDC
// alternate algorithm (using the positive inverse mod R) - see README_REDC.md.
// However, this function does not finalize the output result, and so it does
// not perform the complete REDC by itself.
// Precondition: with theoretical unlimited precision standard multiplication,
//   we require  u < n*R.  The constant R represents the value
//   R = 2^(ma::ma_numeric_limits<T>::digits).
// Returns: the non-finalized montgomery product.
// Implementation: this function is a safer, generic version of the following-
// uint32_t redc_non_finalized(bool& ovf, uint64_t u, uint32_t n, uint32_t invn)
// {
//    using T = uint32_t;
//    using T2 = uint64_t;
//    T m = static_cast<T>(u) * invn;
//    T2 mn = static_cast<T2>(m) * static_cast<T2>(n);
//    T2 t = u - mn;
//    ovf = (u < mn);
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
T REDC_non_finalized2(bool& ovf, T2 u_input, T ulo_input, T n, T inv_n)
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

    // assert(n * inv_n ≡ 1 (mod R))
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<V>(n) * static_cast<V>(inv_n)) == 1);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    V2 u = static_cast<V2>(u_input);
    HPBC_PRECONDITION2(u < static_cast<V2>(n) * R);

    T m = static_cast<T>(static_cast<V>(ulo_input) * static_cast<V>(inv_n));
    V2 mn = static_cast<V2>(m) * static_cast<V2>(n);

    T2 t = static_cast<T2>(u - mn);
    ovf = (u < mn);

    T t_hi = static_cast<T>(static_cast<V2>(t) >> bit_width_T);

    // Postcondition #1
    // ----------------
    // See the proof above Postcondition #1 in RedcLargeR.h's function
    // detail_redc_large::REDC_non_finalized() for proof of this postcondition.
    // For the same reasons given in Postcondition #1 in RedcLargeR.h, we never
    // compute the finalized result, aside from within this postcondition.
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T finalized_result = (ovf) ? static_cast<T>(t_hi + n) : t_hi;
        HPBC_POSTCONDITION2(finalized_result < n);
    }
    // Postcondition #2:  If  n < R/2,  then  0 < t_hi + n < 2*n
    // ---------------------------------------------------------
    // See the proof accompanying Postcondition #2 in RedcLargeR.h's function
    // detail_redc_large::REDC_non_finalized() for proof of this postcondition.
    HPBC_POSTCONDITION2((n < R/2) ? (0 < static_cast<T>(t_hi + n)) &&
                                    (static_cast<T>(t_hi + n) < 2*n) : true);
    // return the non-finalized result
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
  // REDC_non_finalized2().

  static HURCHALLA_FORCE_INLINE
  T REDC(T2 u, T u_lo, T n, T inv_n, FullrangeTag, PrivateAnyTag)
  {
    bool ovf;
    using namespace detail_redc_small;
    T thi = REDC_non_finalized2<T,T2>(ovf, u, u_lo, n, inv_n);
    // We know the following from REDC_non_finalized2()'s Postcondition #1
    T finalized_result = (ovf) ? static_cast<T>(thi + n) : thi;
    HPBC_POSTCONDITION2(finalized_result < n);
    return finalized_result;
  }
  // Since HalfrangeTag inherits from FullrangeTag, it will argument match to
  // the FullrangeTag function above.  I don't see any way to improve upon
  // the FullrangeTag function with a dedicated version for HalfrangeTag.

  static HURCHALLA_FORCE_INLINE
  T REDC(T2 u, T u_lo, T n, T inv_n, QuarterrangeTag, PrivateAnyTag)
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
    T thi = REDC_non_finalized2<T, T2>(ovf, u, u_lo, n, inv_n);
    // Since we have the precondition n < R/4, we know from
    // REDC_non_finalized2's Postcondition #2 that  0 < thi + n < 2*n.
    HPBC_POSTCONDITION2(0 < static_cast<T>(thi + n) &&
                                                 static_cast<T>(thi + n) < 2*n);
    // MontyQuarterRange (and hence QuarterrangeTag) allows any montgomery
    // values that satisfy 0 <= value < 2*n, thus it's safe to return thi + n.
    return static_cast<T>(thi + n);
  }
  // Since SixthrangeTag inherits from QuarterrangeTag, it will argument match
  // to the QuarterrangeTag function above.


  // In principle we could change this file to use the older REDC_non_minimized2
  // (see older git commits of this file) and thereby achieve slightly improved
  // performance for convert_out with FullrangeTag and HalfrangeTag.  See the
  // convert_out discussion in RedcLargeR.h for details.  But since convert_out
  // typically isn't used in performance critical loops, I have elected not to
  // re-introducing/re-enabling the old REDC function.  It seems very unlikely
  // that this function matters enough to be worth that added complexity.
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T inv_n)
  {
    T2 u = static_cast<T2>(x);
    // u = x satisfies REDC_non_finalized2's precondition requiring u < n*R,
    // since u = x < R < n*R.
    bool ovf;
    using namespace detail_redc_small;
    T thi = REDC_non_finalized2(ovf, u, static_cast<T>(u), n, inv_n);
    // REDC_non_finalized2()'s Postcondition #1 guarantees the following
    T finalized_result = (ovf) ? static_cast<T>(thi + n) : thi;
    HPBC_POSTCONDITION2(finalized_result < n);
    return finalized_result;
  }
};


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
