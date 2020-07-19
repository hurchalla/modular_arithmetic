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


#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>


template <typename T>
T brute_modular_pow(T base, T power, T modulus)
{
    namespace ma = hurchalla::modular_arithmetic;
    T result = 1;
    for (T i=0; i<power; ++i)
        result = ma::modular_multiplication_prereduced_inputs(result, base,
                                                                       modulus);
    return result;
}


template <typename T>
void test_modulus(T modulus)
{
    static_cast<void>(modulus);

    namespace ma = hurchalla::modular_arithmetic;

    T base = 0;
    T power = 0;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 0; power = 1;
    EXPECT_TRUE(static_cast<T>(0) == ma::modular_pow(base, power, modulus));
    base = 0; power = 2;
    EXPECT_TRUE(static_cast<T>(0) == ma::modular_pow(base, power, modulus));
    base = 1; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 1; power = 1;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 1; power = 2;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));

    base = static_cast<T>(modulus - 1);
    power = 0;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    power = 1;
    EXPECT_TRUE(base == ma::modular_pow(base, power, modulus));
    power = 2;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    power = 3;
    EXPECT_TRUE(base == ma::modular_pow(base, power, modulus));

    T tmax = ma::ma_numeric_limits<T>::max();
    // make power the largest possible even number
    power = static_cast<T>((tmax/2)*2);
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    --power;  // power should now be odd
    EXPECT_TRUE(base == ma::modular_pow(base, power, modulus));

    base = modulus;
    power = 2;
    EXPECT_TRUE(static_cast<T>(0) == ma::modular_pow(base, power, modulus));
    power = 5;
    EXPECT_TRUE(static_cast<T>(0) == ma::modular_pow(base, power, modulus));

    if (modulus < tmax) {
        base = static_cast<T>(modulus + 1);
        power = 2;
        EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
        power = 5;
        EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    }

    T tmp = static_cast<T>((modulus/4)*4);  // make tmp == 4n for some integer n
    base = static_cast<T>(tmp/2);
    power = 2;
    EXPECT_TRUE(static_cast<T>(0) == ma::modular_pow(base, power, tmp));
}


template <typename T>
void test_modular_pow()
{
    namespace ma = hurchalla::modular_arithmetic;

    // test with a few basic examples first
    T modulus = 13;
    T base = 5;
    T power = 12;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 7; power = 6;
    EXPECT_TRUE(static_cast<T>(12) == ma::modular_pow(base, power, modulus));
    modulus = 14;
    EXPECT_TRUE(static_cast<T>(7) == ma::modular_pow(base, power, modulus));

    base = 5;
    power = 53;
    modulus = 13;
    EXPECT_TRUE(static_cast<T>(5) == ma::modular_pow(base, power, modulus));
    base = 6;
    EXPECT_TRUE(static_cast<T>(2) == ma::modular_pow(base, power, modulus));

    test_modulus(static_cast<T>(13));
    test_modulus(static_cast<T>(14));

    // --------- Test using moduli that are likely edge cases --------

    modulus = 2;
    base = 0; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 0; power = 5;
    EXPECT_TRUE(static_cast<T>(0) == ma::modular_pow(base, power, modulus));
    base = 1; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 31; power = 0;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 1; power = 3;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 17; power = 3;
    EXPECT_TRUE(static_cast<T>(1) == ma::modular_pow(base, power, modulus));
    base = 14; power = 3;
    EXPECT_TRUE(static_cast<T>(0) == ma::modular_pow(base, power, modulus));

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
    TEST(ModularArithmetic, modular_pow) {
        test_modular_pow<uint8_t>();
        test_modular_pow<uint16_t>();
        test_modular_pow<uint32_t>();
        test_modular_pow<uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_pow<__uint128_t>();
#endif

        test_modular_pow<int8_t>();
        test_modular_pow<int16_t>();
        test_modular_pow<int32_t>();
        test_modular_pow<int64_t>();

// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_pow<__int128_t>();
#endif
    }

    TEST(ModularArithmetic, modular_pow_large_exponents) {
        namespace ma = hurchalla::modular_arithmetic;
        // test a couple large exponent cases
        uint32_t base = 81452;
        uint32_t exponent = 113;
        uint32_t modulus = 2951486173u;
        uint32_t result = ma::modular_pow(base, exponent, modulus);
        EXPECT_TRUE(result == brute_modular_pow(base, exponent, modulus));

        base = 81451;
        exponent = 113;
        result = ma::modular_pow(base, exponent, modulus);
        EXPECT_TRUE(result == brute_modular_pow(base, exponent, modulus));

        exponent = 114;
        result = ma::modular_pow(base, exponent, modulus);
        EXPECT_TRUE(result == brute_modular_pow(base, exponent, modulus));
    }
}
