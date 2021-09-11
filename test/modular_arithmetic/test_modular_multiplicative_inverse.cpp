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
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>

// do exhaustive test of all uint8_t?

namespace {


namespace hc = ::hurchalla;

namespace testmmi {
    // For details, see https://en.wikipedia.org/wiki/Greatest_common_divisor
    template <typename T>
    T gcd(T a, T b)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_PRECONDITION(a >= 0);
        HPBC_PRECONDITION(b >= 0);

        while (a != 0) {
            T tmp = a;
            a = static_cast<T>(b % a);
            b = tmp;
        }
        return b;
    }
}  // end namespace testmmi


void exhaustive_test_uint8_t();
void exhaustive_test_uint8_t()
{
    for (std::uint8_t modulus=255; modulus>1; --modulus) {
        for (std::uint8_t a=0; a<modulus; ++a) {
            std::uint8_t inv = hc::modular_multiplicative_inverse(a, modulus);
            if (inv == 0)
                EXPECT_TRUE(1 < testmmi::gcd(a, modulus));
            else
                EXPECT_TRUE(static_cast<std::uint8_t>(1) ==
                    hc::modular_multiplication_prereduced_inputs(a, inv,
                                                                      modulus));
        }
    }
}


template <typename T>
void test_modulus(T modulus)
{
    T a = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = modulus;
    EXPECT_TRUE(static_cast<T>(0) ==
                                hc::modular_multiplicative_inverse(a, modulus));

    T tmax = hc::ut_numeric_limits<T>::max();
    if (modulus < tmax) {
        a = static_cast<T>(modulus + 1);
        EXPECT_TRUE(static_cast<T>(1) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    }

    a = 2;
    T inverse = hc::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < testmmi::gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             hc::modular_multiplication_prereduced_inputs(
                                static_cast<T>(a % modulus), inverse, modulus));

    a = 3;
    inverse = hc::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < testmmi::gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             hc::modular_multiplication_prereduced_inputs(
                                static_cast<T>(a % modulus), inverse, modulus));

    a = static_cast<T>(modulus - 1);
    inverse = hc::modular_multiplicative_inverse(a, modulus);
    EXPECT_TRUE(static_cast<T>(1) ==
             hc::modular_multiplication_prereduced_inputs(a, inverse, modulus));

    a = static_cast<T>(modulus - 2);
    inverse = hc::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < testmmi::gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             hc::modular_multiplication_prereduced_inputs(a, inverse, modulus));

    a = static_cast<T>(modulus/2);
    inverse = hc::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < testmmi::gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             hc::modular_multiplication_prereduced_inputs(a, inverse, modulus));

    a++;
    inverse = hc::modular_multiplicative_inverse(a, modulus);
    if (inverse == 0)
        EXPECT_TRUE(1 < testmmi::gcd(a, modulus));
    else
        EXPECT_TRUE(static_cast<T>(1) ==
             hc::modular_multiplication_prereduced_inputs(a, inverse, modulus));
}



template <typename T>
void test_modular_multiplicative_inverse()
{
    // test with a few basic examples first
    T modulus = 13;
    T a = 5;
    EXPECT_TRUE(static_cast<T>(8) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 7;
    EXPECT_TRUE(static_cast<T>(2) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 4;
    EXPECT_TRUE(static_cast<T>(10) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 17;
    EXPECT_TRUE(static_cast<T>(10) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 14;
    EXPECT_TRUE(static_cast<T>(1) ==
                                hc::modular_multiplicative_inverse(a, modulus));


    // modular_multiplicative_inverse() indicates the inverse doesn't exist by
    // returning 0.  (The inverse doesn't exist when gcd(a, modulus) > 1)
    a = 12;
    modulus = 21;   // a modulus of 21 shares the factor 3 with a.
    EXPECT_TRUE(static_cast<T>(0) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                hc::modular_multiplicative_inverse(a, modulus));

    a = 7;
    modulus = 16;
    EXPECT_TRUE(static_cast<T>(7) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 10;
    EXPECT_TRUE(static_cast<T>(0) ==
                                hc::modular_multiplicative_inverse(a, modulus));

    a = 7;
    modulus = 14;
    EXPECT_TRUE(static_cast<T>(0) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 9;
    EXPECT_TRUE(static_cast<T>(11) ==
                                hc::modular_multiplicative_inverse(a, modulus));

    test_modulus(modulus);
    test_modulus(static_cast<T>(15));

    // --------- Test using moduli that are likely edge cases --------

    modulus = 2;
    a = 0;
    EXPECT_TRUE(static_cast<T>(0) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 1;
    EXPECT_TRUE(static_cast<T>(1) ==
                                hc::modular_multiplicative_inverse(a, modulus));
    a = 5;
    EXPECT_TRUE(static_cast<T>(1) ==
                                hc::modular_multiplicative_inverse(a, modulus));

    modulus = hc::ut_numeric_limits<T>::max();
    test_modulus(modulus);
    modulus--;
    test_modulus(modulus);

    modulus = hc::ut_numeric_limits<T>::max() / 2;
    test_modulus(modulus);
    modulus++;
    test_modulus(modulus);
}



TEST(ModularArithmetic, modular_multiplicative_inverse) {

    exhaustive_test_uint8_t();

    test_modular_multiplicative_inverse<std::uint8_t>();
    test_modular_multiplicative_inverse<std::uint16_t>();
    test_modular_multiplicative_inverse<std::uint32_t>();
    test_modular_multiplicative_inverse<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_modular_multiplicative_inverse<__uint128_t>();
#endif
}


} // end unnamed namespace
