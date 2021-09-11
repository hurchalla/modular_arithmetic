// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_addition.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// Performance note on modular addition vs modular subtraction:
//   You can generally expect that the assembly for a modular subtraction will
// never use more total operations (it will likely use less) than a modular
// addition.  You can also expect that the assembly for modular subtraction will
// never use more total registers (it will likely use less) than modular
// addition.
//   In favor of addition, you can generally expect that the assembly for
// modular addition will never have higher total latency (it will likely be
// lower) than modular subtraction.
//    This comparison is written based on the performance of x86 architecture.
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
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION(modulus>0);
    HPBC_PRECONDITION(a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(b<modulus);   // i.e. the input must be prereduced
    
    T result = detail::impl_modular_addition_prereduced_inputs(a, b, modulus);

    // POSTCONDITION:
    // Returns (a+b)%modulus, performed as if a and b have infinite precision
    // and thus as if (a+b) is never subject to integer overflow.
    HPBC_POSTCONDITION(result<modulus);
    return result;
}


}  // end namespace

#endif
