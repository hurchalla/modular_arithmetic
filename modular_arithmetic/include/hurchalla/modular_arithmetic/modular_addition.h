// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_addition.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// Perfomance notes are given below this function
template <typename T>  HURCHALLA_FORCE_INLINE
T modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(modulus > 0);
    HPBC_PRECONDITION(0<=a && a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(0<=b && b<modulus);   // i.e. the input must be prereduced
    
    T result = detail::impl_modular_addition<T>::call(a, b, modulus);

    // POSTCONDITION:
    // Returns (a+b)%modulus, performed as if a and b have infinite precision
    // and thus as if (a+b) is never subject to integer overflow.
    HPBC_POSTCONDITION(0<=result && result<modulus);
    return result;
}


// --Performance Notes--

// This function can often have lower latency than
// modular_subtraction_prereduced_inputs, and should never have higher latency.
// However, modular_subtraction_prereduced_inputs usually has fewer uops.  If
// you need modular addition and you want to optimize for a low uop count
// (presumably to increase throughput), *and* if you see that 'b' will remain
// constant over many of your modular addition calls (typically due to you
// calling in a loop), then as an option, you could calculate
// negative_b = (modulus - b),  and then, instead of calling
// modular_addition_prereduced_inputs(a, b, modulus),  you can instead call
// modular_subtraction_prereduced_inputs(a, negative_b, modulus).  If you
// calculate negative_b once, and then use it over many calls of
// modular_subtraction_prereduced_inputs, then you will get most of the benefit
// of modular subtraction's lower uop count.

// Performance note for RISC-V (and other uncommon CPU architectures that do not
// have an instruction for conditional move or conditional select):
//   On this architecture, modular addition may perform better when T is signed
// than when it is unsigned.  Specifically, when HURCHALLA_AVOID_CSELECT is
// defined (see hurchalla/util/compiler_macros.h), a signed type may perform
// better; if it is not defined, you should expect no performance difference
// between signed and unsigned.


}  // end namespace

#endif
