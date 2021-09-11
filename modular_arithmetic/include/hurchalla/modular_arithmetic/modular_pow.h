// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/impl_modular_pow.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


template <typename T>
T modular_pow(T base, T exponent, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION(modulus > 1);

    T result = detail::impl_modular_pow(base, exponent, modulus);

    // POSTCONDITION:
    //   Returns the modular exponentiation of base^exponent (mod modulus).
    HPBC_POSTCONDITION(result<modulus);
    return result;
}


}  // end namespace

#endif
