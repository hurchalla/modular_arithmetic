// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// This example is intended for the case that you are not using CMake.
// If you haven't already, you need to follow the steps in the README.md
// for "How to use the library" | "Without CMake"
#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include <iostream>
#include <cassert>
#include <cstdint>

#ifndef NDEBUG
// remove this if you want to allow asserts
// (they're very good for testing and debugging but may drastically slow down
// the library).
#error "Performance warning: asserts are enabled and will slow performance"
#endif

int main()
{
   namespace hc = ::hurchalla;
   int64_t modulus = 333333333;
   int64_t base = 42;
   int64_t exponent = 123456789;

   // ---- Demonstration of modular exponentiation ----

   // Montgomery arithmetic version:
   assert(modulus % 2 == 1);  // montgomery arithmetic always needs odd modulus.
      // first construct a MontgomeryForm object to do Montgomery arithmetic
      // with a particular modulus.
   hc::MontgomeryForm<int64_t> mf(modulus);
      // convert base to its Montgomery representation.
   auto base_montval = mf.convertIn(base);
      // get the pow result in Montgomery representation.
   auto result_montval = mf.pow(base_montval, exponent);
      // convert the Montgomery representation result to normal integer domain.
   int64_t result1 = mf.convertOut(result_montval);


   // Standard arithmetic version:  (note that Montgomery arithmetic is
   // typically faster, and that modular_pow() requires an unsigned type)
   uint64_t result2 = hc::modular_pow(static_cast<uint64_t>(base),
                                      static_cast<uint64_t>(exponent),
                                      static_cast<uint64_t>(modulus));


   std::cout << "Example results for " << base << "^" << exponent
                                           << " (mod " << modulus << ")\n";
   std::cout << "---------\n";
   std::cout << "using Montgomery arithmetic: " << result1 << "\n";
   std::cout << "using standard arithmetic: " << result2 << "\n";

   return 0;
}
