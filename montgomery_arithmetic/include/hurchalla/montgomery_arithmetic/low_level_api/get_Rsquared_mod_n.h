// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_GET_RSQUARED_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_GET_RSQUARED_MOD_N_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic {


// For discussion purposes, let the unlimited precision constant R represent
// R = 2^(ut_numeric_limits<T>::digits).  For example, if T is uint64_t, then
// R = 2^64.

// get_Rsquared_mod_n() computes and returns (R*R) % n.
// You can get the argument inverse_n_modR by calling inverse_mod_r().  You can
// get Rmod_n by calling get_R_mod_n().

template <typename T, class MTAG = FullrangeTag>
T get_Rsquared_mod_n(T n, T inverse_n_modR, T Rmod_n, MTAG = MTAG())
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
    // verify that  n * inverse_n_modR ≡ 1 (mod R)
    HPBC_PRECONDITION2(
       static_cast<T>(static_cast<P>(n) * static_cast<P>(inverse_n_modR)) == 1);

    T rSquaredModN =
             detail::impl_get_Rsquared_mod_n(n, inverse_n_modR, Rmod_n, MTAG());

    HPBC_POSTCONDITION2(rSquaredModN < n);
    return rSquaredModN;
}


}}

#endif
