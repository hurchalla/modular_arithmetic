// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_GET_RSQUARED_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_GET_RSQUARED_MOD_N_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_get_Rsquared_mod_n.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla {


// For discussion purposes, let the unlimited precision constant R represent
// R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t, then
// R = 1<<64.

// get_Rsquared_mod_n() computes and returns (R*R) % n.
// You can get the argument inverse_n_modR by calling inverse_mod_r().  You can
// get Rmod_n by calling get_R_mod_n().

template <typename T>
T get_Rsquared_mod_n(T n, T inverse_n_modR, T Rmod_n)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(n % 2 == 1);  // REDC requires an odd modulus.
    HPBC_PRECONDITION2(n > 1);
    using P = typename safely_promote_unsigned<T>::type;
    // verify that  n * inverse_n_modR ≡ 1 (mod R)
    HPBC_PRECONDITION2(
       static_cast<T>(static_cast<P>(n) * static_cast<P>(inverse_n_modR)) == 1);

    T rSquaredModN =
               detail::impl_get_Rsquared_mod_n::call(n, inverse_n_modR, Rmod_n);

    HPBC_POSTCONDITION2(rSquaredModN < n);
    return rSquaredModN;
}


}

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
