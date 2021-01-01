// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_LLAPI_REDC_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_LLAPI_REDC_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_REDC.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic {


// This REDC function provides the alternate REDC algorithm described at
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

// For discussion purposes below, let the unlimited precision constant R
// represent  R = 2^(ut::ut_numeric_limits<T>::digits).  For example, if T is
// uint64_t, then R = 2^64.

// When calling REDC, you must either accept the default for the template
// parameter MTAG or use one of the structs from monty_tag_structs.h.  Likewise
// for PTAG, you must must either accept the default or use one of the structs
// from optimization_tag_structs.h.
// For details on the benefits and requirements of different MTAGs and PTAGs,
// see respectively, monty_tag_structs.h and optimization_tag_structs.h.

template <typename T, class MTAG = FullrangeTag, class PTAG = LowlatencyTag>
HURCHALLA_FORCE_INLINE
T REDC(T u_hi, T u_lo, T n, T inv_n, MTAG = MTAG(), PTAG = PTAG())
{
    namespace ut = hurchalla::util;
    static_assert(ut::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut::ut_numeric_limits<T>::is_signed), "");
    static_assert(ut::ut_numeric_limits<T>::is_modulo, "");
    if (HPBC_PRECONDITION2_MACRO_IS_ACTIVE) {
        // Using HalfrangeTag requires that the modulus n < R/2.
        // QuarterrangeTag requires n < R/4.  SixthrangeTag requires n < R/6.
        T Rdiv2 = static_cast<T>(
                   static_cast<T>(1) << (ut::ut_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2((std::is_same<MTAG, HalfrangeTag>::value) ? 
                           n < Rdiv2 : true);
        HPBC_PRECONDITION2((std::is_same<MTAG, QuarterrangeTag>::value) ? 
                           n < (Rdiv2/2) : true);
        HPBC_PRECONDITION2((std::is_same<MTAG, SixthrangeTag>::value) ? 
                           n < (Rdiv2/3) : true);
    }
    HPBC_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_PRECONDITION2(n > 1);
    using P = typename ut::safely_promote_unsigned<T>::type;
    // verify that  n * inv_n ≡ 1 (mod R)
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<P>(n) * static_cast<P>(inv_n)) == 1);
    HPBC_PRECONDITION2(u_hi < n);  // verify that (u_hi*R + u_lo) < n*R

    T result = detail::impl_REDC(u_hi, u_lo, n, inv_n, MTAG(), PTAG());

    HPBC_POSTCONDITION2(
        (std::is_same<MTAG, FullrangeTag>::value ||
         std::is_same<MTAG, HalfrangeTag>::value) ?
        result < n :
        0 < result && result < 2*n);
    return result;
}


// This overload of REDC provides a parameter resultIsZero that lets you know
// if the return value is congruent to zero (mod n).  REDC() can determine this
// more efficiently than would be possible for a user, for certain MTAG types.
template <typename T, class MTAG = FullrangeTag, class PTAG = LowlatencyTag>
HURCHALLA_FORCE_INLINE
T REDC(T u_hi, T u_lo, T n, T inv_n, bool& resultIsZero, MTAG = MTAG(),
                                                                  PTAG = PTAG())
{
    // All preconditions and postconditions are the same as REDC() above, since
    // we just call it from here.
    T result = REDC(u_hi, u_lo, n, inv_n, MTAG(), PTAG());
    resultIsZero = detail::isZeroRedcResult(result, n, MTAG());
    return result;
}


}} // end namespace

#endif