
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_POW_H__INCLUDED


#include "hurchalla/modular_arithmetic/standard/internal/impl_modular_pow.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


// Interface/contract.
template <typename T>
T modular_pow(T base, T exponent, T modulus)  noexcept
{
    static_assert(std::is_unsigned<T>::value, "");  //T unsigned integral type
    precondition(modulus > 1);
    // Postcondition:
    //   Returns the modular exponentiation of base^exponent (mod modulus).

    return impl_modular_pow(base, exponent, modulus);
}


}}  // end namespace

#endif
