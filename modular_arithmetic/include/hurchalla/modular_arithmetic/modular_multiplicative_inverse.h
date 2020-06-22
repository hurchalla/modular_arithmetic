// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATIVE_INVERSE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATIVE_INVERSE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/impl_modular_multiplicative_inverse.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace modular_arithmetic {


// Note: Calling with a < modulus slightly improves performance.
// [The multiplicative inverse is an integer > 0 and < modulus, such that
//    a * multiplicative_inverse == 1 (mod modulus).   It is a unique number,
//    but it exists if and only if 'a' and 'modulus' are coprime.]
template <typename T>
T modular_multiplicative_inverse(T a, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(a >= 0);
    HPBC_PRECONDITION(modulus > 1);

    T inverse = impl_modular_multiplicative_inverse(a, modulus);

    HPBC_POSTCONDITION(inverse >= 0 && inverse < modulus);
    //POSTCONDITION: Returns 0 if the inverse does not exist. Otherwise returns
    //   the value of the inverse (which is never 0, given that modulus>1).
    HPBC_POSTCONDITION(inverse == 0 || modular_multiplication_prereduced_inputs(
                                           a % modulus, inverse, modulus) == 1);
    return inverse;
}


}}	// end namespace

#endif
