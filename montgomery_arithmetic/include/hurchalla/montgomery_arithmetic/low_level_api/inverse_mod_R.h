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
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

namespace hurchalla {


// For discussion purposes, let type UP be a conceptually unlimited precision
// unsigned integer type, and let the unlimited precision constant R represent
// R = (UP)1 << ut_numeric_limits<T>::digits.  Equivalently,
// R = (UP)ut_numeric_limits<T>::max + 1.  For example, if T is uint64_t, we
// would have R = (UP)1 << 64.

// Returns the integer x satisfying  x*a ≡ 1 (mod R)
// This function is constexpr when compiling for std=c++14 or higher
template <typename T>
HURCHALLA_CPP14_CONSTEXPR
T inverse_mod_R(T a)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    HPBC_CLOCKWORK_CONSTEXPR_PRECONDITION(a % 2 == 1);

    T inv = detail::impl_inverse_mod_R::call<T,ut_numeric_limits<T>::digits>(a);

    // guarantee inv*a ≡ 1 (mod R)
    using P = typename safely_promote_unsigned<T>::type;
    HPBC_CLOCKWORK_CONSTEXPR_POSTCONDITION(static_cast<T>(1) ==
                       static_cast<T>(static_cast<P>(inv) * static_cast<P>(a)));
    return inv;
}


} // end namespace

#endif
