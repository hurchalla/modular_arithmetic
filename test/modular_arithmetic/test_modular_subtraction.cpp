// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


// We'll define HURCHALLA_ALLOW_INLINE_ASM_MODSUB here in order to make modular
// subtraction use an inline asm function version if it is available.
// Internally, this inline asm function will also call the generic template
// function version of modular subtraction inside a postcondition, in order to
// make sure that the asm result is correct.  Of course postcondition checks
// must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to undefine NDEBUG, which is why we undef
// NDEBUG here too.
#undef HURCHALLA_ALLOW_INLINE_ASM_MODSUB
#define HURCHALLA_ALLOW_INLINE_ASM_MODSUB 1
#undef NDEBUG


#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>


template <typename T>
void test_modulus(T modulus)
{
    namespace ma = hurchalla::modular_arithmetic;

    T a = 0;
    T b = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    a = 0; b = 1;
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));
    a = 1; b = 1;
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));

    a = 0; b = static_cast<T>(modulus - 1);
    EXPECT_TRUE(static_cast<T>(1) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(b, b, modulus));

    a = 1; b = static_cast<T>(modulus - 1);
    EXPECT_TRUE(static_cast<T>(2) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(modulus - 2) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));

    a = 0; b = static_cast<T>(modulus - 2);
    EXPECT_TRUE(static_cast<T>(2) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(modulus - 2) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));

    a = static_cast<T>(modulus/2);
    b = static_cast<T>(a + 1);
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(a, a, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(b, b, modulus));

    b++;
    EXPECT_TRUE(static_cast<T>(modulus - 2) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(2) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));
    a++;
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));

    a = static_cast<T>(modulus/2 - 1);
    b = static_cast<T>(a + 1);
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(a, a, modulus));
}


template <typename T>
void test_modular_subtraction()
{
    namespace ma = hurchalla::modular_arithmetic;

    // test with a few basic examples first
    T modulus = 13;
    T a = 5;
    T b = 12;
    EXPECT_TRUE(static_cast<T>(6) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(7) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(b, b, modulus));
    a = 7; b = 6;
    EXPECT_TRUE(static_cast<T>(1) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(12) ==
                      ma::modular_subtraction_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(b, b, modulus));

    test_modulus(modulus);
    test_modulus(static_cast<T>(14));

    // --------- Test using moduli that are likely edge cases --------

    modulus = 1;
    a = 0; b = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                      ma::modular_subtraction_prereduced_inputs(a, b, modulus));

    modulus = ma::ma_numeric_limits<T>::max();
    test_modulus(modulus);
    modulus--;
    test_modulus(modulus);

    modulus = ma::ma_numeric_limits<T>::max() / 2;
    test_modulus(modulus);
    modulus++;
    test_modulus(modulus);
}



namespace {
    TEST(ModularArithmetic, modular_subtraction) {
        test_modular_subtraction<uint8_t>();
        test_modular_subtraction<uint16_t>();
        test_modular_subtraction<uint32_t>();
        test_modular_subtraction<uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_subtraction<__uint128_t>();
#endif

        test_modular_subtraction<int8_t>();
        test_modular_subtraction<int16_t>();
        test_modular_subtraction<int32_t>();
        test_modular_subtraction<int64_t>();

// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_subtraction<__int128_t>();
#endif
    }
}
