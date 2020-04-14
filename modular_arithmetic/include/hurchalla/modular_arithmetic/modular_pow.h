
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H__INCLUDED


#include "hurchalla/modular_arithmetic/internal/impl_modular_pow.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace modular_arithmetic {


// Interface/contract.
template <typename T>
T modular_pow(T base, T exponent, T modulus)
{
    static_assert(std::numeric_limits<T>::is_integer &&
                 !(std::numeric_limits<T>::is_signed), "");
    precondition(modulus > 1);
    // Postcondition:
    //   Returns the modular exponentiation of base^exponent (mod modulus).

    return impl_modular_pow(base, exponent, modulus);
}


}}  // end namespace

#endif
