// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_GET_R_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_GET_R_MOD_N_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// For discussion purposes, let the unlimited precision constant R represent
// R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t, then
// R = 1<<64.

// Compute R % n
template <typename T>
T get_R_mod_n(T n)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    // Assign a tmp T variable rather than directly using the intermediate
    // expression, in order to avoid a negative value (and a wrong answer)
    // in cases where 'n' would be promoted to type 'int'.
    T tmp = static_cast<T>(static_cast<T>(0) - n);
    // Compute R%n.  For example, if R==1<<64, arithmetic wraparound behavior
    // of the unsigned integral type T results in (0 - n) representing
    // ((1<<64)-n).  Thus, rModN = R%n == (1<<64)%n == ((1<<64)-n)%n == (0-n)%n
    T rModN = static_cast<T>(tmp % n);
    // Since n is odd and > 1, n does not divide R==1<<x.  Thus, rModN != 0.

    HPBC_POSTCONDITION2(0 < rModN && rModN < n);
    return rModN;
}


}

#endif
