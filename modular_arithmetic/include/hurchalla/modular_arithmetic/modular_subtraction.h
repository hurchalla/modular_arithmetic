// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_SUBTRACTION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_SUBTRACTION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace modular_arithmetic {


template <typename T>
T modular_subtraction_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(modulus>0);
    HPBC_PRECONDITION(a>=0 && a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(b>=0 && b<modulus);   // i.e. the input must be prereduced

    // POSTCONDITION:
    // Returns (a-b)%modulus, performed as if a and b are infinite precision
    // signed ints and thus as if (a-b) is never subject to wraparound/overflow.

    // We want essentially-  result = (a-b < 0) ? a-b+modulus : a-b
    //    But for unsigned type T, (a-b < 0) is always false.  So instead we use
    T tmp = static_cast<T>(a-b);
    T result = (a<b) ? static_cast<T>(modulus+tmp) : tmp;
    return result;
}


}}  // end namespace

#endif
