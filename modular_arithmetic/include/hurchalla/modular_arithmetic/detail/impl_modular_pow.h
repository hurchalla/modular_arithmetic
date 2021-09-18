// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
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
#if 0
    T result = (exponent & 1) ? base : 1;
#else
    T result = 1;
    HURCHALLA_CMOV(exponent & 1, result, base);
#endif
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
