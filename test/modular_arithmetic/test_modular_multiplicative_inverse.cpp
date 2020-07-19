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


#include "hurchalla/modular_arithmetic/modular_multiplicative_inverse.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>


// do exhaustive test of all uint8_t?


// For details, see https://en.wikipedia.org/wiki/Greatest_common_divisor
template <typename T>
T gcd(T a, T b)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(a >= 0);
    HPBC_PRECONDITION(b >= 0);

    while (a != 0) {
        T tmp = a;
        a = static_cast<T>(b % a);
        b = tmp;
    }
    return b;
}


void exhaustive_test_uint8_t();
void exhaustive_test_uint8_t()
{
    namespace ma = hurchalla::modular_arithmetic;
    for (uint8_t modulus = 255; modulus>1; --modulus) {
        for (uint8_t a=0; a<modulus; ++a) {
            uint8_t inverse = ma::modular_multiplicative_inverse(a, modulus);
            if (inverse == 0)
                EXPECT_TRUE(1 < gcd(a, modulus));
            else
                EXPECT_TRUE(static_cast<uint8_t>(1) ==
                    ma::modular_multiplication_prereduced_inputs(a, inverse,
                                                                      modulus));
        }
    }
}


template <typename T>
void test_modulus(T modulus)
{
    namespace ma = hurchalla::modular_arithmetic;

    T a = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = modulus;
    EXPECT_TRUE(static_cast<T>(0) ==
                                ma::modular_multiplicative_inverse(a, modulus));

    T tmax = ma::ma_numeric_limits<T>::max();
    if (modulus < tmax) {
        a = static_cast<T>(modulus + 1);
        EXPECT_TRUE(static_cast<T>(1) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    }

    a = 2;
    T inverse = ma::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             ma::modular_multiplication_prereduced_inputs(
                                static_cast<T>(a % modulus), inverse, modulus));

    a = 3;
    inverse = ma::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             ma::modular_multiplication_prereduced_inputs(
                                static_cast<T>(a % modulus), inverse, modulus));

    a = static_cast<T>(modulus - 1);
    inverse = ma::modular_multiplicative_inverse(a, modulus);
    EXPECT_TRUE(static_cast<T>(1) ==
             ma::modular_multiplication_prereduced_inputs(a, inverse, modulus));

    a = static_cast<T>(modulus - 2);
    inverse = ma::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             ma::modular_multiplication_prereduced_inputs(a, inverse, modulus));

    a = static_cast<T>(modulus/2);
    inverse = ma::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             ma::modular_multiplication_prereduced_inputs(a, inverse, modulus));

    a++;
    inverse = ma::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             ma::modular_multiplication_prereduced_inputs(a, inverse, modulus));
}



template <typename T>
void test_modular_multiplicative_inverse()
{
    namespace ma = hurchalla::modular_arithmetic;

    // test with a few basic examples first
    T modulus = 13;
    T a = 5;
    EXPECT_TRUE(static_cast<T>(8) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 7;
    EXPECT_TRUE(static_cast<T>(2) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 4;
    EXPECT_TRUE(static_cast<T>(10) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 17;
    EXPECT_TRUE(static_cast<T>(10) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 14;
    EXPECT_TRUE(static_cast<T>(1) ==
                                ma::modular_multiplicative_inverse(a, modulus));


    // modular_multiplicative_inverse() indicates the inverse doesn't exist by
    // returning 0.  (The inverse doesn't exist when gcd(a, modulus) > 1)
    a = 12;
    modulus = 21;   // a modulus of 21 shares the factor 3 with a.
    EXPECT_TRUE(static_cast<T>(0) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                ma::modular_multiplicative_inverse(a, modulus));

    a = 7;
    modulus = 16;
    EXPECT_TRUE(static_cast<T>(7) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 10;
    EXPECT_TRUE(static_cast<T>(0) ==
                                ma::modular_multiplicative_inverse(a, modulus));

    a = 7;
    modulus = 14;
    EXPECT_TRUE(static_cast<T>(0) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 9;
    EXPECT_TRUE(static_cast<T>(11) ==
                                ma::modular_multiplicative_inverse(a, modulus));

    test_modulus(modulus);
    test_modulus(static_cast<T>(15));

    // --------- Test using moduli that are likely edge cases --------

    modulus = 2;
    a = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                ma::modular_multiplicative_inverse(a, modulus));
    a = 5;
    EXPECT_TRUE(static_cast<T>(1) ==
                                ma::modular_multiplicative_inverse(a, modulus));

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
    TEST(ModularArithmetic, modular_multiplicative_inverse) {

        exhaustive_test_uint8_t();

        test_modular_multiplicative_inverse<uint8_t>();
        test_modular_multiplicative_inverse<uint16_t>();
        test_modular_multiplicative_inverse<uint32_t>();
        test_modular_multiplicative_inverse<uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_multiplicative_inverse<__uint128_t>();
#endif

        test_modular_multiplicative_inverse<int8_t>();
        test_modular_multiplicative_inverse<int16_t>();
        test_modular_multiplicative_inverse<int32_t>();
        test_modular_multiplicative_inverse<int64_t>();

// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_modular_multiplicative_inverse<__int128_t>();
#endif
    }
}
