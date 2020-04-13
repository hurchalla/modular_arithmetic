
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H__INCLUDED


#include "hurchalla/modular_arithmetic/internal/impl_modular_multiplication.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


// Interface/contract.
template <typename T>
T modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(std::is_unsigned<T>::value, "");  //T unsigned integral type
    precondition(modulus>0);
    precondition(a<modulus);    // i.e. the input must be prereduced
    precondition(b<modulus);    // i.e. the input must be prereduced
    // Postcondition: Returns (a*b)%modulus, theoretically calculated at
    //                infinite precision to avoid overflow.

    return impl_modular_multiplication_prereduced_inputs(a, b, modulus);
}


}} // end namespace

#endif
