// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// Returns the modular exponentiation of base^exponent (mod modulus).
// For details, see http://en.wikipedia.org/wiki/Modular_exponentiation
// note: uses a static member function to disallow ADL.
struct impl_modular_pow {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T base, T exponent, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus > 1);

    if (base >= modulus)
       base = static_cast<T>(base % modulus);
/*
    // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
    // Algorithm 14.76, original unoptimized version
    T result = 1;
    while (exponent > 0)
    {
       if (exponent & 1)
          result= modular_multiplication_prereduced_inputs(result,base,modulus);
       exponent = exponent >> 1;
       base = modular_multiplication_prereduced_inputs(base, base, modulus);
    }
*/
    // slightly optimized version
       // result = (exponent & 1) ? base : 1;
    T result = conditional_select((exponent & 1), base, static_cast<T>(1));
    while (exponent > 1)
    {
       exponent = static_cast<T>(exponent >> 1);
       base = modular_multiplication_prereduced_inputs(base, base, modulus);
       if (exponent & 1)
          result= modular_multiplication_prereduced_inputs(result,base,modulus);
    }
    return static_cast<T>(result);
  }
};


}}  // end namespace

#endif  // include guard
