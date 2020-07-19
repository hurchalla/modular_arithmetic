// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/modular_multiplicative_inverse.h"
#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/montgomery_arithmetic/detail/MontgomeryValue.h"

#include "gtest/gtest.h"
#include <cstdint>



struct hardcoded_test_67 {
    template <typename T, typename M>
    void test() const
    {
        T modulus = 67;
        // this test isn't applicable if the modulus is larger than the max
        // modulus that the template param M (the MontyType) allows.
        if (modulus > M::max_modulus())
            return;

        M mf(modulus);

        using V = typename decltype(mf)::V;
        V x = mf.convertIn(60);
        V y = mf.convertIn(13);

        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 6);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 20);
        EXPECT_TRUE(mf.getCanonicalForm(mf.add(x,y)) ==
                                          mf.getCanonicalForm(mf.convertIn(6)));
        EXPECT_TRUE(mf.getUnityValue() == mf.convertIn(1));
        EXPECT_TRUE(mf.getZeroValue() == mf.convertIn(0));
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                                     mf.convertIn(static_cast<T>(modulus - 1)));
        EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == 43);
        EXPECT_TRUE(mf.convertOut(mf.square(y)) == 35);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 1)) == 13);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 2)) == 35);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 5)) == 46);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 7)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 8)) == 26);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 11)) == 38);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 12)) == 25);
    }
};



#ifdef __SIZEOF_INT128__
struct hardcoded_test_uint128 {
    template <typename T, typename M>
    void test() const
    {
        // For this test, we require the template param T == __uint128_t.
        // Using a template param T is a bit of a hack for compatibility
        // with test_default() and test_explicit()
        static_assert(std::is_same<T, __uint128_t>::value, "");

        T modulus =
                 hurchalla::modular_arithmetic::ma_numeric_limits<T>::max() - 2;
        // this test isn't applicable if the modulus is larger than the max
        // modulus that the template param M (the MontyType) allows.
        if (modulus > M::max_modulus())
            return;

        M mf(modulus);
    
        using V = typename decltype(mf)::V;
        V x = mf.convertIn(modulus - 1);
        V y = mf.convertIn(2);
    
        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 3);
        EXPECT_TRUE(mf.getCanonicalForm(mf.add(x,y)) ==
                                        mf.getCanonicalForm(mf.convertIn(1)));
        EXPECT_TRUE(mf.getUnityValue() == mf.convertIn(1));
        EXPECT_TRUE(mf.getZeroValue() == mf.convertIn(0));
        EXPECT_TRUE(mf.getNegativeOneValue() == mf.convertIn(modulus - 1));
        EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == modulus - 2);
        EXPECT_TRUE(mf.convertOut(mf.square(x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 1)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 2)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 10)) == 1024);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 128)) == 3);
    }
};
#endif



template <typename T, typename F>
void test_default()
{
    namespace ma = hurchalla::montgomery_arithmetic;
    F f{};
    f.template test<T, ma::MontgomeryForm<T>>();
}

template <typename T, typename F>
void test_explicit()
{
    namespace ma = hurchalla::montgomery_arithmetic;
    F f{};
    f.template test<T, ma::MontgomeryForm<T,ma::MontyWrappedStandardMath<T>>>();
    f.template test<T, ma::MontgomeryForm<T,ma::MontyFullRange<T>>>();
    f.template test<T, ma::MontgomeryForm<T,ma::MontyHalfRange<T>>>();
    f.template test<T, ma::MontgomeryForm<T,ma::MontyQuarterRange<T>>>();
    f.template test<T, ma::MontgomeryForm<T,ma::MontySqrtRange<T>>>();
}




namespace {
    // Test basic montgomery functions with a hardcoded modulus=67
    // and hardcoded arguments for the functions.  We'll test every
    // permutation of types possible without using extended
    // precision libraries:

    TEST(MontgomeryArithmetic67Default, int8_t) {
        test_default<int8_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Default, int16_t) {
        test_default<int16_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Default, int32_t) {
        test_default<int32_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Default, int64_t) {
        test_default<int64_t, hardcoded_test_67>();
    }
#ifdef __SIZEOF_INT128__
    TEST(MontgomeryArithmetic67Default, __int128_t) {
        test_default<__int128_t, hardcoded_test_67>();
    }
#endif

    TEST(MontgomeryArithmetic67Default, uint8_t) {
        test_default<uint8_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Default, uint16_t) {
        test_default<uint16_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Default, uint32_t) {
        test_default<uint32_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Default, uint64_t) {
        test_default<uint64_t, hardcoded_test_67>();
    }
#ifdef __SIZEOF_INT128__
    TEST(MontgomeryArithmetic67Default, __uint128_t) {
        test_default<__uint128_t, hardcoded_test_67>();
    }
#endif

    TEST(MontgomeryArithmetic67Explicit, uint8_t) {
        test_explicit<uint8_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Explicit, uint16_t) {
        test_explicit<uint16_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Explicit, uint32_t) {
        test_explicit<uint32_t, hardcoded_test_67>();
    }
    TEST(MontgomeryArithmetic67Explicit, uint64_t) {
        test_explicit<uint64_t, hardcoded_test_67>();
    }
#ifdef __SIZEOF_INT128__
    TEST(MontgomeryArithmetic67Explicit, __uint128_t) {
        test_explicit<__uint128_t, hardcoded_test_67>();
    }
#endif

#ifdef __SIZEOF_INT128__
    TEST(MontgomeryArithmeticFull, uint128_tests) {
        test_default<__uint128_t, hardcoded_test_uint128>();
        test_explicit<__uint128_t, hardcoded_test_uint128>();
    }
#endif
}
