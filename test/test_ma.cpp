
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


template <typename T>
T brute_modular_pow(T a, T power, T modulus)
{
    namespace ma = hurchalla::modular_arithmetic;
    T result = 1;
    for (T i=0; i<power; ++i)
        result = ma::modular_multiplication_prereduced_inputs(result, a,
                                                              modulus);
    return result;
}



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
    TEST(ModularArithmetic, basic_examples) {
        using T = uint32_t;
        namespace ma = hurchalla::modular_arithmetic;
        T a = 5;
        T b = 12;
        T modulus = 13;

        T result = ma::modular_multiplication_prereduced_inputs(a, b, modulus);
        EXPECT_TRUE(result == static_cast<T>(8));
        result = ma::modular_addition_prereduced_inputs(a, b, modulus);
        EXPECT_TRUE(result == static_cast<T>(4));
        result = ma::modular_subtraction_prereduced_inputs(a, b, modulus);
        EXPECT_TRUE(result == static_cast<T>(6));

        T base = 5;
        T exponent = 53;
        result = ma::modular_pow(base, exponent, modulus);
        EXPECT_TRUE(result == static_cast<T>(5));
        base = 81452;
        exponent = 4013;
        modulus = 2951486173u;
        result = ma::modular_pow(base, exponent, modulus);
        EXPECT_TRUE(result == brute_modular_pow(base, exponent, modulus));

        modulus = 13;
        result = ma::modular_multiplicative_inverse(a, modulus);
        EXPECT_TRUE(result == static_cast<T>(8));
        // The inverse doesn't exist if gcd(a, modulus) > 1.
        // A return value of 0 indicates the inverse doesn't exist.
        modulus = 21;   // shares the factor 3 with b
        result = ma::modular_multiplicative_inverse(b, modulus);
        EXPECT_TRUE(result == static_cast<T>(0));
    }


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
