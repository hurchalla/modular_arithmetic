// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// Alternatively, please consider using the montgomery_multiplication class
// MontgomeryForm (specifically its multiply function) instead of this function
// modular_multiplication_prereduced_inputs().  If you are heavily using modular
// multiplication in your code, there's a decent chance that montgomery
// multiplication will improve performance- often significantly.  It always
// requires an odd modulus though.

template <typename T>
T modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION(modulus>0);
    HPBC_PRECONDITION(a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(b<modulus);   // i.e. the input must be prereduced

    T result =
           detail::impl_modular_multiplication_prereduced_inputs(a, b, modulus);

    // POSTCONDITION: Returns (a*b)%modulus, theoretically calculated at
    //                infinite precision to avoid overflow.
    HPBC_POSTCONDITION(result<modulus);
    return result;
}


} // end namespace

#endif
