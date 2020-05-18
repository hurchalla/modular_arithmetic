
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED


#include "hurchalla/modular_arithmetic/internal/impl_modular_multiplication.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace modular_arithmetic {


// Interface/contract.
template <typename T>
T modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(modulus>0);
    HPBC_PRECONDITION(a>=0 && a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(b>=0 && b<modulus);   // i.e. the input must be prereduced

    // POSTCONDITION: Returns (a*b)%modulus, theoretically calculated at
    //                infinite precision to avoid overflow.

    return impl_modular_multiplication_prereduced_inputs(a, b, modulus);
}


}} // end namespace

#endif
