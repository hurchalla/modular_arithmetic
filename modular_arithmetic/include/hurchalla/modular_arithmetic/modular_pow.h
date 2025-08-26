// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/impl_modular_pow.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

namespace hurchalla {


// Alternatively, please consider using the MontgomeryForm class member function
// pow() instead of this function modular_pow().  There's an excellent chance
// that you will achieve much better perfomance using MontgomeryForm's pow -
// though note that MontgomeryForm can only be used if your modulus is odd.

template <typename T, typename U>
T modular_pow(T base, U exponent, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");
    HPBC_CLOCKWORK_API_PRECONDITION(modulus > 1);

    T result = detail::impl_modular_pow::call(base, exponent, modulus);

    // POSTCONDITION:
    //  Returns the modular exponentiation of base to the exponent (mod modulus)
    HPBC_CLOCKWORK_POSTCONDITION(result<modulus);
    return result;
}


}  // end namespace

#endif
