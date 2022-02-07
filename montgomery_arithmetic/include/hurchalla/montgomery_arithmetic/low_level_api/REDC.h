// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_LLAPI_REDC_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_LLAPI_REDC_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/ImplRedc.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla {


// This file is the API for the alternate REDC algorithm described at
// https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/README_REDC.md
// This alternate version of the REDC algorithm differs in small but important
// ways from Peter Montgomery's original 1985 paper "Modular multiplication
// without trial division".  From the point of view of a caller, the most
// important distinction is that this version requires the positive inverse for
// one of its arguments rather than the negative inverse (which was required by
// the original/traditional REDC algorithm).  We provide the alternate version
// instead of the traditional version, because it improves efficiency both in
// terms of latency and number of instructions.  See README_REDC.md for the
// details.
//
// For discussion purposes below, let the unlimited precision constant R
// represent  R = 2^(ut::ut_numeric_limits<T>::digits).  For example, if T is
// uint64_t, then R = 2^64.


//   REDC_standard() returns the standard and normally expected value from REDC,
// which is the least residue (modulo the modulus).  In other words,
// REDC_standard() guarantees that 0 <= return_value < modulus.
//   When calling REDC_standard(), you must either accept the default for the
// template parameter PTAG or use one of the structs from
// optimization_tag_structs.h  - see optimization_tag_structs.h for the
// benefits and requirements of different PTAGs.
template <typename T, class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
T REDC_standard(T u_hi, T u_lo, T n, T inv_n, PTAG = PTAG())
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");

    HPBC_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_PRECONDITION2(n > 1);
    using P = typename safely_promote_unsigned<T>::type;
    // verify that  n * inv_n ≡ 1 (mod R)
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_PRECONDITION2(u_hi < n);  // verify that (u_hi*R + u_lo) < n*R

    T result = detail::RedcStandard<T>::call(u_hi, u_lo, n, inv_n, PTAG());

    HPBC_POSTCONDITION2(result < n);
    return result;
}


// REDC_incomplete() is "incomplete" in that the return value is not finalized
// to the least residue mod n (i.e. it might not satisfy 0 <= return_value < n;
// see the postconditions).  The reason you might wish to use this function is
// that it can provide better performance than the standard REDC in some
// situations (examples further below).
// The isNegative parameter allows you to effectively interpret the return value
// as a signed two's complement value, with a bit width one larger than T.  For
// example, if type T is uint64_t, then you can view the return value as beng a
// 65 bit signed (two's complement) value, in which the (implied) 65th bit is
// represented by isNegative.  This signed result is guaranteed to be congruent
// (mod n) to the result from REDC_standard(); more specifically, if isNegative
// is true, this function's return value equals  n + REDC_standard(), and if
// isNegative is false, the return value equals REDC_standard().
// ---
// As an example of an optimization that REDC_incomplete() allows, we can use it
// to optimize Montgomery multiplication when the modulus 'n' is less than R/4.
// For such a case, we can unconditionally add 'n' to the return value from
// REDC_incomplete(), and we can then use this sum "as-is" as an input to a
// Montgomery multiplication.  There is no need to minimize this input value to
// the least residue (i.e. to  0 <= sum < n), because Montgomery multiplication
// with inputs x and y only requires that  u = x*y < n*R,  and we will always be
// able to satisfy this requirement when using the described sum:
// For any values 'a' and 'b' returned by REDC_incomplete() when n < R/4,  let
// x=a+n  and  y=b+n.  We know by our REDC_incomplete() postcondition that if
// n < R/2 then  0 < a+n < 2*n  and  0 < b+n < 2*n.  Thus the multiplication
// u = x*y == (a+n)*(b+n) < (2*n)*(2*n) == (4*n)*n < (4*n)*(R/4) == n*R,  which
// means u < n*R, as required.  For more details, see also section 5 from
// "Montgomery's Multiplication Technique: How to Make It Smaller and Faster"
// https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps
// You can see also MontyQuarterRange.h in this library, which is a class
// that requires n < R/4 and is optimized in this way.
// Another example of a different optimization enabled by REDC_incomplete is
// MontyHalfrange.h.
template <typename T> HURCHALLA_FORCE_INLINE
T REDC_incomplete(bool& isNegative, T u_hi, T u_lo, T n, T inv_n)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");

    HPBC_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_PRECONDITION2(n > 1);
    using P = typename safely_promote_unsigned<T>::type;
    // verify that  n * inv_n ≡ 1 (mod R)
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_PRECONDITION2(u_hi < n);  // verify that (u_hi*R + u_lo) < n*R

    T result = detail::RedcIncomplete::call(isNegative, u_hi, u_lo, n, inv_n);

    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T finalized_result = (isNegative) ? static_cast<T>(result + n) : result;
        HPBC_POSTCONDITION2(finalized_result ==
                                           REDC_standard(u_hi, u_lo, n, inv_n));
        HPBC_POSTCONDITION2(0 <= finalized_result && finalized_result < n);
    }
    // If  n < R/2,  then  0 < result + n < 2*n
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T Rdiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
        HPBC_POSTCONDITION2((n < Rdiv2) ? (0 < static_cast<T>(result + n)) &&
                                (static_cast<T>(result + n) < 2*n) : true);
    }
    return result;
}


} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
