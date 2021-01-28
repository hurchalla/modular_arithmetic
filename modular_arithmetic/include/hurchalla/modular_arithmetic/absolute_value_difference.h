// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


template <typename T>
T absolute_value_difference(T a, T b)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(a >= 0);
    HPBC_PRECONDITION(b >= 0);

    T result = detail::impl_absolute_value_difference(a, b);

    // POSTCONDITION:
    // This function returns absolute_value(a-b).  A simple potential
    // implementation could be  result = (a>b) ? a-b : b-a.
    HPBC_POSTCONDITION(result<=a || result<=b);
    return result;
}


}  // end namespace

#endif
