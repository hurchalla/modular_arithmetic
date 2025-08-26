// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// This example is intended for the case that you are using CMake.
// If you haven't already, you need to follow the steps in the README.md
// for "How to use the library" | "With CMake"
#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include <iostream>
#include <cassert>
#include <cstdint>


int main()
{
   namespace hc = ::hurchalla;

   // you could use any integer type that the compiler supports
   // (including __uint128_t)
   using T = uint64_t;

   T modulus = 333333333;
   T base = 42;
   T exponent = 123456789;

   // ---- Demonstration of modular exponentiation ----

   // Montgomery arithmetic version:
   assert(modulus % 2 == 1);  // montgomery arithmetic always needs odd modulus.
      // First construct a MontgomeryForm object to do Montgomery arithmetic
      // with the modulus we chose.
   hc::MontgomeryForm<T> mf(modulus);
      // Convert base to its Montgomery representation.
   auto mont_base = mf.convertIn(base);
      // Get the pow result in Montgomery representation.
   auto mont_result = mf.pow(mont_base, exponent);
      // Convert the Montgomery representation result to normal integer domain.
   T result1 = mf.convertOut(mont_result);


   // Standard arithmetic version:  (note that Montgomery arithmetic is
   // usually much faster)
   T result2 = hc::modular_pow(base, exponent, modulus);


   std::cout << "Example results for " << base << "^" << exponent
                                           << " (mod " << modulus << ")\n";
   std::cout << "---------\n";
   std::cout << "using Montgomery arithmetic: " << result1 << "\n";
   std::cout << "using standard arithmetic: " << result2 << "\n";

   return 0;
}
