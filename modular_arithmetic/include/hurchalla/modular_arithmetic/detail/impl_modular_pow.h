// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace modular_arithmetic {


// Code review/testing notes: Analytically, modular_pow appears perfect
//    (assuming modular_multiplication_prereduced_inputs() is correct).
//    Empirically, impossible to exhaustively test, but passed all tests.
// Adapted from pseudocode @ http://en.wikipedia.org/wiki/Modular_exponentiation
template <typename T>
T impl_modular_pow(T base, T exponent, T modulus)
{
   static_assert(ma_numeric_limits<T>::is_integer, "");
   static_assert(!(ma_numeric_limits<T>::is_signed), "");
   HPBC_PRECONDITION2(modulus > 1);
   if (base >= modulus)
      base = base % modulus;
/*
   // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
   // Algorithm 14.76, original unoptimized version
   T result = 1;
   while (exponent > 0)
   {
      if (exponent & 1)
         result = modular_multiplication_prereduced_inputs(result,base,modulus);
      exponent = exponent >> 1;
      base = modular_multiplication_prereduced_inputs(base, base, modulus);
   }
*/
   // slightly optimized version
   T result;
   if (exponent & 1)
      result = base;
   else
      result = 1;
   while (exponent > 1)
   {
      exponent = exponent >> 1;
      base = modular_multiplication_prereduced_inputs(base, base, modulus);
      if (exponent & 1)
         result = modular_multiplication_prereduced_inputs(result,base,modulus);
   }
   return result;
}


}}  // end namespace

#endif  // include guard
