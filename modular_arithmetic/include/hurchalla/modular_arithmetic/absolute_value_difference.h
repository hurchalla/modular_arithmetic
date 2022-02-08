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


template <typename T>
T absolute_value_difference(T a, T b)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    T result = detail::impl_absolute_value_difference<T>::call(a, b);

    HPBC_POSTCONDITION(result == ((a>b) ? a-b : b-a));
    HPBC_POSTCONDITION(result<=a || result<=b);
    return result;
}


}  // end namespace

#endif
