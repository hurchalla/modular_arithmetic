// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// We'll undefine HURCHALLA_DISALLOW_INLINE_ASM_MODMUL here in order to make
// modular multiplication use an inline asm function version if it is available.
// This shouldn't be strictly necessary, since there's no reason this macro
// would be defined at this point, and by default modular multiplication uses
// inline asm (if available) unless this macro is defined.
// Internally, the inline asm function will also call the generic template
// function version of modular multiplication inside a postcondition, in order
// to make sure that the asm result is correct.  Of course postcondition checks
// must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to undefine NDEBUG, which is why we undef
// NDEBUG here too.
#undef HURCHALLA_DISALLOW_INLINE_ASM_MODMUL
#undef NDEBUG


#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>

namespace {


namespace hc = ::hurchalla;

template <typename T>
T brute_modular_pow(T base, T power, T modulus)
{
    T result = 1;
    for (T i=0; i<power; ++i)
        result = hc::modular_multiplication_prereduced_inputs(result, base,
                                                                       modulus);
    return result;
}


template <typename T>
void test_modulus(T modulus)
{
    static_cast<void>(modulus);

    T base = 0;
    T power = 0;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 0; power = 1;
    EXPECT_TRUE(static_cast<T>(0) == hc::modular_pow(base, power, modulus));
    base = 0; power = 2;
    EXPECT_TRUE(static_cast<T>(0) == hc::modular_pow(base, power, modulus));
    base = 1; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 1; power = 1;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 1; power = 2;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));

    base = static_cast<T>(modulus - 1);
    power = 0;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    power = 1;
    EXPECT_TRUE(base == hc::modular_pow(base, power, modulus));
    power = 2;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    power = 3;
    EXPECT_TRUE(base == hc::modular_pow(base, power, modulus));

    T tmax = hc::ut_numeric_limits<T>::max();
    // make power the largest possible even number
    power = static_cast<T>((tmax/2)*2);
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    --power;  // power should now be odd
    EXPECT_TRUE(base == hc::modular_pow(base, power, modulus));

    base = modulus;
    power = 2;
    EXPECT_TRUE(static_cast<T>(0) == hc::modular_pow(base, power, modulus));
    power = 5;
    EXPECT_TRUE(static_cast<T>(0) == hc::modular_pow(base, power, modulus));

    if (modulus < tmax) {
        base = static_cast<T>(modulus + 1);
        power = 2;
        EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
        power = 5;
        EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    }

    T tmp = static_cast<T>((modulus/4)*4);  // make tmp == 4n for some integer n
    base = static_cast<T>(tmp/2);
    power = 2;
    EXPECT_TRUE(static_cast<T>(0) == hc::modular_pow(base, power, tmp));
}


template <typename T>
void test_modular_pow()
{
    // test with a few basic examples first
    T modulus = 13;
    T base = 5;
    T power = 12;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 7; power = 6;
    EXPECT_TRUE(static_cast<T>(12) == hc::modular_pow(base, power, modulus));
    modulus = 14;
    EXPECT_TRUE(static_cast<T>(7) == hc::modular_pow(base, power, modulus));

    base = 5;
    power = 53;
    modulus = 13;
    EXPECT_TRUE(static_cast<T>(5) == hc::modular_pow(base, power, modulus));
    base = 6;
    EXPECT_TRUE(static_cast<T>(2) == hc::modular_pow(base, power, modulus));

    test_modulus(static_cast<T>(13));
    test_modulus(static_cast<T>(14));

    // --------- Test using moduli that are likely edge cases --------

    modulus = 2;
    base = 0; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 0; power = 5;
    EXPECT_TRUE(static_cast<T>(0) == hc::modular_pow(base, power, modulus));
    base = 1; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 31; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 1; power = 3;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 17; power = 3;
    EXPECT_TRUE(static_cast<T>(1) == hc::modular_pow(base, power, modulus));
    base = 14; power = 3;
    EXPECT_TRUE(static_cast<T>(0) == hc::modular_pow(base, power, modulus));

    modulus = hc::ut_numeric_limits<T>::max();
    test_modulus(modulus);
    modulus--;
    test_modulus(modulus);

    modulus = hc::ut_numeric_limits<T>::max() / 2;
    test_modulus(modulus);
    modulus++;
    test_modulus(modulus);
}



TEST(ModularArithmetic, modular_pow) {
    test_modular_pow<std::uint8_t>();
    test_modular_pow<std::uint16_t>();
    test_modular_pow<std::uint32_t>();
    test_modular_pow<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_modular_pow<__uint128_t>();
#endif
}

TEST(ModularArithmetic, modular_pow_large_exponents) {
    // test a couple large exponent cases
    std::uint32_t base = 81452;
    std::uint32_t exponent = 113;
    std::uint32_t modulus = 2951486173u;
    std::uint32_t result = hc::modular_pow(base, exponent, modulus);
    EXPECT_TRUE(result == brute_modular_pow(base, exponent, modulus));

    base = 81451;
    exponent = 113;
    result = hc::modular_pow(base, exponent, modulus);
    EXPECT_TRUE(result == brute_modular_pow(base, exponent, modulus));

    exponent = 114;
    result = hc::modular_pow(base, exponent, modulus);
    EXPECT_TRUE(result == brute_modular_pow(base, exponent, modulus));
}


} // end unnamed namespace
