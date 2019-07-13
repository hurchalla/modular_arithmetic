
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATIVE_INVERSE_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATIVE_INVERSE_H__INCLUDED


#include "hurchalla/modular_arithmetic/standard/internal/impl_modular_multiplicative_inverse.h"
#include "hurchalla/modular_arithmetic/standard/modular_multiplication.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


// Note: Calling with a < modulus slightly improves performance.
// [The multiplicative inverse is an integer > 0 and < modulus, such that
//    a * multiplicative_inverse == 1 (mod modulus).   It is a unique number,
//    but it exists if and only if 'a' and 'modulus' are coprime.]
template <typename T>
T modular_multiplicative_inverse(T a, T modulus)  noexcept
{
    static_assert(std::is_unsigned<T>::value, "");  //T unsigned integral type
    precondition(modulus>1);

    T inverse = impl_modular_multiplicative_inverse(a, modulus);

    postcondition(inverse < modulus);
    //Postcondition: Returns 0 if the inverse does not exist. Otherwise returns
    //   the value of the inverse (which is never 0, given that modulus>1).
    postcondition(inverse == 0 || modular_multiplication_prereduced_inputs(
                                           a % modulus, inverse, modulus) == 1);
    return inverse;
}


}}	// end namespace

#endif
