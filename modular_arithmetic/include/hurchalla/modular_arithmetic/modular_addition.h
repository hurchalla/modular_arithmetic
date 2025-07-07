// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/impl_modular_addition.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla {


// Perfomance recommendations are given below this function
template <typename T, class PTAG = LowlatencyTag>  HURCHALLA_FORCE_INLINE
T modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(std::is_same<PTAG,LowlatencyTag>::value ||
                  std::is_same<PTAG,LowuopsTag>::value, "");
    HPBC_PRECONDITION(modulus > 0);
    HPBC_PRECONDITION(0<=a && a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(0<=b && b<modulus);   // i.e. the input must be prereduced
    
    T result = detail::impl_modular_addition<T,PTAG>::call(a, b, modulus);

    // POSTCONDITION:
    // Returns (a+b)%modulus, performed as if a and b have infinite precision
    // and thus as if (a+b) is never subject to integer overflow.
    HPBC_POSTCONDITION(0<=result && result<modulus);
    return result;
}


// --Performance Suggestions--

// Note on the optional PTAG (performance tag) template parameter:
//    This parameter does not affect functionality in any way; it affects only
// performance.
//    If you specify PTAG, you must choose either LowuopsTag or LowlatencyTag.
//    If you want low latency, you should usually choose LowlatencyTag, provided
// that you can see that neither 'modulus' nor 'b' were recently set or modified
// prior to your call - note that a "recent" modify could be on a prior loop
// iteration.  If either of them was recently modified, the compiler will
// usually be unable to provide any low latency benefit over LowuopsTag.  Note
// that LowlatencyTag will typically use more uops and create more pressure on
// the ALU than LowuopsTag.
//    You should usually prefer LowuopsTag if you want to minimize the uop count
// and ALU pressure (presumably for higher throughput), or if you see that
// either 'b' or 'modulus' was set/modified close to your call of this function
// (note that "close" could be on a prior loop iteration).  LowuopsTag generally
// provides a lower uop count and lower ALU pressure than LowlatencyTag.  We can
// note though that it is possible for LowlatencyTag to effectively have similar
// uop count and ALU pressure as LowuopsTag, if the compiler can loop hoist
// LowlatencyTag's extra instruction(s) involving 'b' and 'modulus'.

// Performance note for RISC-V (and other uncommon CPU architectures that do not
// have an instruction for conditional move or conditional select):
//   On this architecture, modular addition may perform better when T is signed
// than when it is unsigned.  Specifically, when HURCHALLA_AVOID_CSELECT is
// defined (see hurchalla/util/compiler_macros.h), a signed type may perform
// better; if it is not defined, you should expect no performance difference
// between signed and unsigned.


}  // end namespace

#endif
