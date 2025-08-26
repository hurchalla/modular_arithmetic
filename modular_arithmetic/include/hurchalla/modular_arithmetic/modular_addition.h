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
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

namespace hurchalla {


// Perfomance notes are given below this function
template <typename T>  HURCHALLA_FORCE_INLINE
T modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    HPBC_CLOCKWORK_API_PRECONDITION(modulus > 0);
    HPBC_CLOCKWORK_API_PRECONDITION(0<=a && a<modulus);   // i.e. the input must be prereduced
    HPBC_CLOCKWORK_API_PRECONDITION(0<=b && b<modulus);   // i.e. the input must be prereduced
    
    T result = detail::impl_modular_addition<T>::call(a, b, modulus);

    // POSTCONDITION:
    // Returns (a+b)%modulus, performed as if a and b have infinite precision
    // and thus as if (a+b) is never subject to integer overflow.
    HPBC_CLOCKWORK_POSTCONDITION(0<=result && result<modulus);
    return result;
}


// --Performance Notes--

// For this function to be able to complete at its lowest latency, you will need
// to ensure in your calling code (if possible) that neither 'b' nor 'modulus'
// was recently changed (or set) prior to your call of this function - note that
// "recently" could be on a prior loop iteration.  Generally speaking, if either
// 'b' or 'modulus' was changed on an immediately preceding (modular) arithmetic
// instruction, or if one of those two variables was otherwise changed
// immediately beforehand, then usually this function will need one more cycle
// to complete than it would need at its ideal lowest latency.
// If you wish to maximize throughput rather than minimize latency, then you may
// find modular subtraction to be useful - modular subtraction by default has
// fewer uops than modular addition (note that subtraction never has lower
// latency).  Given that modular subtraction by default uses fewer uops, if you
// need to do modular additions and you want to optimize for a low uop count,
// *and* if you see that 'b' will remain constant over many of your modular
// addition calls (typically due to you calling in a loop), then as an option,
// you could calculate:
// negative_b = (b == 0) ? 0 : (modulus - b),  and then, instead of calling
// modular_addition_prereduced_inputs(a, b, modulus),  you can instead call
// modular_subtraction_prereduced_inputs(a, negative_b, modulus).  If you
// calculate negative_b once, and then use it over many calls of
// modular_subtraction_prereduced_inputs, then potentially modular subtraction's
// lowered uop count might increase your overall throughput slightly.

// Performance note for RISC-V (and other uncommon CPU architectures that do not
// have an instruction for conditional move or conditional select):
//   On this architecture, modular addition may perform better when T is signed
// than when it is unsigned.  Specifically, when HURCHALLA_AVOID_CSELECT is
// defined (see hurchalla/util/compiler_macros.h), a signed type may perform
// better; if it is not defined, you should expect no performance difference
// between signed and unsigned.


}  // end namespace

#endif
