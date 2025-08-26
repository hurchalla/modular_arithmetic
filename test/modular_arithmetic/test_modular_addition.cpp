// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// Strictly for testing purposes, we'll define HURCHALLA_ALLOW_INLINE_ASM_ALL
// here in order to make modular addition use an inline asm function version if
// it is available.  Internally, this inline asm function will also call the
// generic template function version of modular addition inside a postcondition,
// in order to make sure that the asm result is correct. Of course postcondition
// checks must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to define HURCHALLA_CLOCKWORK_ENABLE_ASSERTS,
// which is why we do so here.  This is all strictly for testing purposes.
#undef HURCHALLA_ALLOW_INLINE_ASM_ALL
#define HURCHALLA_ALLOW_INLINE_ASM_ALL 1

#ifndef HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#  define HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#endif


#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>

namespace {


namespace hc = ::hurchalla;

template <typename T>
void test_modulus(T modulus)
{
    EXPECT_TRUE(modulus > 2); // if this fails, this test file has a bug

    T a = 0;
    T b = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    a = 0; b = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));
    a = 1; b = 1;
    EXPECT_TRUE(static_cast<T>(2) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));

    a = 0; b = static_cast<T>(modulus - 1);
    EXPECT_TRUE(b == hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(b == hc::modular_addition_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(modulus - 2) ==
                         hc::modular_addition_prereduced_inputs(b, b, modulus));

    a = 1; b = static_cast<T>(modulus - 1);
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));

    a = static_cast<T>(modulus/2);
    b = static_cast<T>(modulus - a);
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));

    b++;
    EXPECT_TRUE(static_cast<T>(1) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));
    a++;
    EXPECT_TRUE(static_cast<T>(2) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(2) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));

    a = static_cast<T>(modulus/2 - 1);
    b = static_cast<T>(modulus - a - 2);
    EXPECT_TRUE(static_cast<T>(modulus - 2) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(modulus - 2) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));
    a++;
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));
}


template <typename T>
void test_modular_addition()
{
    // test with a few basic examples first
    T modulus = 13;
    T a = 5;
    T b = 12;
    EXPECT_TRUE(static_cast<T>(4) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(4) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(11) ==
                         hc::modular_addition_prereduced_inputs(b, b, modulus));
    a = 7; b = 6;
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(12) ==
                         hc::modular_addition_prereduced_inputs(b, b, modulus));

    test_modulus(modulus);
    test_modulus(static_cast<T>(14));

    // --------- Test using moduli that are likely edge cases --------

    modulus = 1;
    a = 0; b = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                         hc::modular_addition_prereduced_inputs(a, b, modulus));

    modulus = hc::ut_numeric_limits<T>::max();
    test_modulus(modulus);
    modulus--;
    test_modulus(modulus);

    modulus = hc::ut_numeric_limits<T>::max() / 2;
    test_modulus(modulus);
    modulus++;
    test_modulus(modulus);
}



TEST(ModularArithmetic, modular_addition) {
    test_modular_addition<std::uint8_t>();
    test_modular_addition<std::uint16_t>();
    test_modular_addition<std::uint32_t>();
    test_modular_addition<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_modular_addition<__uint128_t>();
#endif

    test_modular_addition<std::int8_t>();
    test_modular_addition<std::int16_t>();
    test_modular_addition<std::int32_t>();
    test_modular_addition<std::int64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_modular_addition<__int128_t>();
#endif
}


} // end unnamed namespace
