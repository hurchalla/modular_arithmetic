// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


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


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
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
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    a = 0; b = 1;
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));
    a = 1; b = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));

    a = 2; b = 3;
    EXPECT_TRUE(static_cast<T>(6) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(6) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(4) ==
                   ma::modular_multiplication_prereduced_inputs(a, a, modulus));

    a = 0; b = static_cast<T>(modulus - 1);
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                   ma::modular_multiplication_prereduced_inputs(b, b, modulus));

    a = 1; b = static_cast<T>(modulus - 1);
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(modulus - 1) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));

    a = static_cast<T>(modulus - 1);
    b = static_cast<T>(modulus - 2);
    EXPECT_TRUE(static_cast<T>(2) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(2) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));

    a = static_cast<T>(modulus - 2);
    b = static_cast<T>(modulus - 3);
    EXPECT_TRUE(static_cast<T>(6) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(6) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));

    T tmp = static_cast<T>((modulus/4)*4);  // make tmp == 4n for some integer n
    a = static_cast<T>(tmp/2);
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(a, a, tmp));

    tmp = static_cast<T>((modulus/2)*2);
    a = static_cast<T>(tmp/2);
    b = static_cast<T>(6);
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, tmp));
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, tmp));

    b = static_cast<T>(5);
    EXPECT_TRUE(a == ma::modular_multiplication_prereduced_inputs(a, b, tmp));
    EXPECT_TRUE(a == ma::modular_multiplication_prereduced_inputs(b, a, tmp));
}


template <typename T>
void test_modular_multiplication()
{
    namespace ma = hurchalla::modular_arithmetic;

    // test with a few basic examples first
    T modulus = 13;
    T a = 5;
    T b = 12;
    EXPECT_TRUE(static_cast<T>(8) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(8) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));
    EXPECT_TRUE(static_cast<T>(12) ==
                   ma::modular_multiplication_prereduced_inputs(a, a, modulus));
    EXPECT_TRUE(static_cast<T>(1) ==
                   ma::modular_multiplication_prereduced_inputs(b, b, modulus));

    modulus = 14;
    a = 7;
    b = 8;
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(b, a, modulus));

    test_modulus(modulus);
    test_modulus(static_cast<T>(15));

    // --------- Test using moduli that are likely edge cases --------

    modulus = 1;
    a = 0; b = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                   ma::modular_multiplication_prereduced_inputs(a, b, modulus));

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
    TEST(ModularArithmetic, modular_multiplication) {
        test_modular_multiplication<uint8_t>();
        test_modular_multiplication<uint16_t>();
        test_modular_multiplication<uint32_t>();
        test_modular_multiplication<uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_multiplication<__uint128_t>();
#endif

        test_modular_multiplication<int8_t>();
        test_modular_multiplication<int16_t>();
        test_modular_multiplication<int32_t>();
        test_modular_multiplication<int64_t>();

// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_multiplication<__int128_t>();
#endif
    }
}
