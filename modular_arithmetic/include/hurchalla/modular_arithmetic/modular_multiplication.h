// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_MULTIPLICATION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace modular_arithmetic {


// Interface/contract.
template <typename T>
T modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    namespace ut = hurchalla::util;
    static_assert(ut::ut_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(modulus>0);
    HPBC_PRECONDITION(a>=0 && a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(b>=0 && b<modulus);   // i.e. the input must be prereduced

    T result =
           detail::impl_modular_multiplication_prereduced_inputs(a, b, modulus);

    // POSTCONDITION: Returns (a*b)%modulus, theoretically calculated at
    //                infinite precision to avoid overflow.
    HPBC_POSTCONDITION(0<=result && result<modulus);
    return result;
}


}} // end namespace

#endif
