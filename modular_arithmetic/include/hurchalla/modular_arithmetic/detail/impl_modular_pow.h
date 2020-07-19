// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/modular_arithmetic/traits/extensible_make_unsigned.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace modular_arithmetic {


// See http://en.wikipedia.org/wiki/Modular_exponentiation
// Code review notes: Analytically, modular_pow() appears correct (assuming that
//     modular_multiplication_prereduced_inputs() is correct).
template <typename T>
HURCHALLA_FORCE_INLINE 
T impl_modular_pow(T base_t, T exponent_t, T modulus_t)
{
   static_assert(ma_numeric_limits<T>::is_integer, "");
   HPBC_PRECONDITION2(modulus_t > 1);
   HPBC_PRECONDITION2(base_t >= 0);
   HPBC_PRECONDITION2(exponent_t >= 0);

   using U = typename extensible_make_unsigned<T>::type;
   U base = static_cast<U>(base_t);
   U exponent = static_cast<U>(exponent_t);
   U modulus = static_cast<U>(modulus_t);

   if (base >= modulus)
      base = static_cast<U>(base % modulus);
/*
   // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
   // Algorithm 14.76, original unoptimized version
   U result = 1;
   while (exponent > 0)
   {
      if (exponent & 1)
         result = modular_multiplication_prereduced_inputs(result,base,modulus);
      exponent = exponent >> 1;
      base = modular_multiplication_prereduced_inputs(base, base, modulus);
   }
*/
   // slightly optimized version
   U result;
   if (exponent & 1)
      result = base;
   else
      result = 1;
   while (exponent > 1)
   {
      exponent = static_cast<U>(exponent >> 1);
      base = modular_multiplication_prereduced_inputs(base, base, modulus);
      if (exponent & 1)
         result = modular_multiplication_prereduced_inputs(result,base,modulus);
   }
   return static_cast<T>(result);
}


}}  // end namespace

#endif  // include guard
