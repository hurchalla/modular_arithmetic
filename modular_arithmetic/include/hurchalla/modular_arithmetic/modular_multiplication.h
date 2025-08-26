// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

namespace hurchalla {


// Alternatively, please consider using the montgomery multiplication class
// MontgomeryForm (specifically its multiply function) instead of this function
// modular_multiplication_prereduced_inputs().  If you are heavily using modular
// multiplication in your code, there's a very good chance that montgomery
// multiplication will improve performance- often significantly.  It always
// requires an odd modulus though.

template <typename T>
T modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_CLOCKWORK_API_PRECONDITION(modulus>0);
    HPBC_CLOCKWORK_API_PRECONDITION(a<modulus);   // i.e. the input must be prereduced
    HPBC_CLOCKWORK_API_PRECONDITION(b<modulus);   // i.e. the input must be prereduced

    T result = detail::impl_modular_multiplication<T>::call(a, b, modulus);

    // POSTCONDITION: Returns (a*b)%modulus, theoretically calculated at
    //                infinite precision to avoid overflow.
    HPBC_CLOCKWORK_POSTCONDITION(result<modulus);
    return result;
}

// You may find the function modular_multiplication_has_slow_perf() to be useful
// when you have a calculation that seems borderline as to whether standard
// modular multiplication or montgomery multiplication would perform better, in
// general across systems.  You can use this function to help you choose at
// compile-time whether you will perform a montgomery or standard modmult (e.g.
// by calling this function in a constexpr if).
template <typename T>
HURCHALLA_FORCE_INLINE constexpr bool modular_multiplication_has_slow_perf()
{
    return detail::impl_modular_multiplication<T>::has_slow_perf();
}


} // end namespace

#endif
