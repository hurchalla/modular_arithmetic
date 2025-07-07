// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// This is an exhaustive test of modular addition using modular subtraction to
// verify the addition results.  The test uses only type uint8_t, in order to
// make it computationaly feasible.


// Strictly for testing purposes, we'll define HURCHALLA_ALLOW_INLINE_ASM_ALL
// here in order to make modular addition use an inline asm function version if
// it is available.  Internally, this inline asm function will also call the
// generic template function version of modular addition inside a postcondition,
// in order to make sure that the asm result is correct.  Of course
// postcondition checks must be enabled for this check to occur - the easiest
// way to ensure postconditions are enabled is to undefine NDEBUG, which is why
// we undef NDEBUG here too.  The same applies to modular subtraction as well.
#undef HURCHALLA_ALLOW_INLINE_ASM_ALL
#define HURCHALLA_ALLOW_INLINE_ASM_ALL 1
#undef NDEBUG


#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "gtest/gtest.h"
#include <cstdint>

namespace {


TEST(ModularArithmetic, modular_addition_with_subtraction) {
    namespace hc = ::hurchalla;
    using T = std::uint8_t;

    for (T modulus=255; modulus>0; --modulus) {
        for (T a=0; a<modulus; ++a) {
            for (T b=0; b<modulus; ++b) {
                T sum = hc::modular_addition_prereduced_inputs(a, b, modulus);
                EXPECT_TRUE(a == (hc::modular_subtraction_prereduced_inputs<T,hc::LowlatencyTag>(sum, b, modulus)));
                EXPECT_TRUE(a == (hc::modular_subtraction_prereduced_inputs<T,hc::LowuopsTag>(sum, b, modulus)));
            }
        }
    }
}


} // end unnamed namespace
