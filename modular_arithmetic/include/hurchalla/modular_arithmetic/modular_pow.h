// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/impl_modular_pow.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace modular_arithmetic {


// Interface/contract.
template <typename T>
T modular_pow(T base, T exponent, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(modulus > 1);
    HPBC_PRECONDITION(base >= 0);
    HPBC_PRECONDITION(exponent >= 0);

    T result = impl_modular_pow(base, exponent, modulus);

    // POSTCONDITION:
    //   Returns the modular exponentiation of base^exponent (mod modulus).
    HPBC_POSTCONDITION(0<=result && result<modulus);
    return result;
}


}}  // end namespace

#endif
