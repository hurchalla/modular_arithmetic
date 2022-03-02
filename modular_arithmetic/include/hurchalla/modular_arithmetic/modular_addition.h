// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_addition.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// Performance note #1 (signed vs. unsigned T):
//   On some systems, this function may perform better when T is signed than
// when it is unsigned.  Specifically, when HURCHALLA_AVOID_CSELECT is defined
// (see hurchalla/util/compiler_macros.h) a signed type can perform better; if
// it is not defined you should expect no performance difference between signed
// and unsigned.
//
// Performance note #2 (modular addition vs modular subtraction):
//   You can generally expect that the assembly for a modular subtraction will
// never use more total operations (it will likely use less) than a modular
// addition.  You can also expect that the assembly for modular subtraction will
// never use more total registers (it will likely use less) than modular
// addition.
//   In favor of addition, you can generally expect that the assembly for
// modular addition will never have higher total latency (it will likely be
// lower) than modular subtraction.
//   This comparison is written based on the performance of x86 architecture.
// If possible, these relative performance characteristic will be maintained
// across architectures and into the future, but these characteristics are not
// guaranteed.
//
// If you will be adding or subtracting the same value over every iteration of
// a loop, then if you wish, you can convert a loop's add into a subtract, or a
// subtract into an add.  You can do this by negating the relevant operand prior
// to loop start.  For example, if you have the loop:
//      for (int i=0; i<total; ++i)
//          x = modular_addition_prereduced_inputs(x, b, modulus);
// you can change it to:
//      T c = modular_subtraction_prereduced_inputs(0, b, modulus);  // negate b
//      for (int i=0; i<total; ++i)
//          x = modular_subtraction_prereduced_inputs(x, c, modulus);


template <typename T>
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


}  // end namespace

#endif
