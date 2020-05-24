
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/impl_modular_pow.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace modular_arithmetic {


// Interface/contract.
template <typename T>
T modular_pow(T base, T exponent, T modulus)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION(modulus > 1);

    // POSTCONDITION:
    //   Returns the modular exponentiation of base^exponent (mod modulus).

    return impl_modular_pow(base, exponent, modulus);
}


}}  // end namespace

#endif
