// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATIVE_INVERSE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATIVE_INVERSE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/impl_modular_multiplicative_inverse.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// Note: Calling with a < modulus slightly improves performance.
// [The multiplicative inverse is an integer > 0 and < modulus, such that
//    a * multiplicative_inverse == 1 (mod modulus).   It is a unique number,
//    but it exists if and only if 'a' and 'modulus' are coprime.]
template <typename T>
T modular_multiplicative_inverse(T a, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION(modulus > 1);

    T inverse = detail::impl_modular_multiplicative_inverse::call(a, modulus);

    HPBC_POSTCONDITION(inverse < modulus);
    //POSTCONDITION: Returns 0 if the inverse does not exist. Otherwise returns
    //   the value of the inverse (which is never 0, given that modulus>1).
    HPBC_POSTCONDITION(inverse == 0 ||
                       ::hurchalla::modular_multiplication_prereduced_inputs(
                           static_cast<T>(a % modulus), inverse, modulus) == 1);
    return inverse;
}


}	// end namespace

#endif
