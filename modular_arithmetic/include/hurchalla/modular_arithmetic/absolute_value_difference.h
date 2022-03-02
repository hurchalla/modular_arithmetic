// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// Performance note:
//   On some systems, this function may perform better when T is signed than
// when it is unsigned.  Specifically, when HURCHALLA_AVOID_CSELECT is defined
// (see hurchalla/util/compiler_macros.h) a signed type can perform better; if
// it is not defined you should expect no performance difference between signed
// and unsigned.


template <typename T>
T absolute_value_difference(T a, T b)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(a >= 0);
    HPBC_PRECONDITION(b >= 0);

    T result = detail::impl_absolute_value_difference<T>::call(a, b);

    HPBC_POSTCONDITION(result >= 0);
    HPBC_POSTCONDITION(result == ((a>b) ? a-b : b-a));
    return result;
}


}  // end namespace

#endif
