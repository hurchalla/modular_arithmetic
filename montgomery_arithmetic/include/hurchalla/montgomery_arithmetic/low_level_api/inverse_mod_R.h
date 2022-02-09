// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_INVERSE_MOD_R_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_INVERSE_MOD_R_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/detail/impl_inverse_mod_R.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// For discussion purposes, let the unlimited precision constant R equal
// 1<<(ut_numeric_limits<T>::digits). For example when T is uint64_t, R = 1<<64.

// Returns the integer x satisfying  x*a ≡ 1 (mod R)
// This function is constexpr when compiling for std=c++14 or higher
template <typename T>
HURCHALLA_CPP14_CONSTEXPR
T inverse_mod_R(T a)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    HPBC_CONSTEXPR_PRECONDITION(a % 2 == 1);

    T inv = detail::impl_inverse_mod_R::call<T,ut_numeric_limits<T>::digits>(a);

    // guarantee inv*a ≡ 1 (mod R)
    using P = typename safely_promote_unsigned<T>::type;
    HPBC_CONSTEXPR_POSTCONDITION(static_cast<T>(1) ==
                       static_cast<T>(static_cast<P>(inv) * static_cast<P>(a)));
    return inv;
}


} // end namespace

#endif
