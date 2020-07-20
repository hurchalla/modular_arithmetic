// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


// For more complete testing we'll define the following macros,
//    HURCHALLA_ALLOW_INLINE_ASM_MODADD
//    HURCHALLA_ALLOW_INLINE_ASM_MODSUB
//    HURCHALLA_ALLOW_INLINE_ASM_MONTMUL
// which will get MontgomeryForm to use any helper inline asm functions that
// are available.  Internally, these inline asm functions will also call their
// corresponding generic template helper functions inside a postcondition, in
// order to make sure that the asm result is correct.  Of course postcondition
// checks must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to undefine NDEBUG, which is why we undef
// NDEBUG here too.
#undef HURCHALLA_ALLOW_INLINE_ASM_MODADD
#define HURCHALLA_ALLOW_INLINE_ASM_MODADD 1
#undef HURCHALLA_ALLOW_INLINE_ASM_MODSUB
#define HURCHALLA_ALLOW_INLINE_ASM_MODSUB 1
#undef HURCHALLA_ALLOW_INLINE_ASM_MONTMUL
#define HURCHALLA_ALLOW_INLINE_ASM_MONTMUL 1
#undef NDEBUG


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "gtest/gtest.h"
#include <cstdint>


template <typename M>
void test_mf_general_checks(M& mf, typename M::T_type a, typename M::T_type b)
{
    namespace ma = hurchalla::modular_arithmetic;

    using T = typename M::T_type;
    using V = typename M::V;
    T modulus = mf.getModulus();
    V x = mf.convertIn(a);
    V y = mf.convertIn(b);

    T reference_sum = ma::modular_addition_prereduced_inputs(a,b,modulus);
    EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == reference_sum);
    EXPECT_TRUE(mf.getCanonicalForm(mf.add(x,y)) ==
                              mf.getCanonicalForm(mf.convertIn(reference_sum)));

    EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) ==
                        ma::modular_subtraction_prereduced_inputs(b,a,modulus));
    EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) ==
                        ma::modular_subtraction_prereduced_inputs(a,b,modulus));

    EXPECT_TRUE(mf.getUnityValue() == mf.getCanonicalForm(mf.convertIn(1)));
    EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalForm(mf.convertIn(0)));
    EXPECT_TRUE(mf.getNegativeOneValue() ==
              mf.getCanonicalForm(mf.convertIn(static_cast<T>(modulus-1))));

    T ref_product = ma::modular_multiplication_prereduced_inputs(a,b,modulus);
    EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == ref_product);
    EXPECT_TRUE(mf.convertOut(mf.multiply(y,x)) == ref_product);

    EXPECT_TRUE(mf.convertOut(mf.square(x)) ==
                     ma::modular_multiplication_prereduced_inputs(a,a,modulus));
    EXPECT_TRUE(mf.convertOut(mf.square(y)) ==
                     ma::modular_multiplication_prereduced_inputs(b,b,modulus));

    EXPECT_TRUE(mf.convertOut(mf.pow(y,0)) == 1);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,1)) == b);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,4)) == ma::modular_pow<T>(b,4,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,13)) ==
                                              ma::modular_pow<T>(b,13,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,17)) ==
                                              ma::modular_pow<T>(b,17,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,127)) ==
                                             ma::modular_pow<T>(b,127,modulus));
}


template <typename M>
void test_MontgomeryForm()
{
    using T = typename M::T_type;
    using V = typename M::V;

    // Try a basic test case first that is valid for all possible Monty types,
    // even M == MontySqrtRange<uint8_t>.
    {
        T modulus = 13;
        M mf(modulus);
        V x = mf.convertIn(6);
        V y = mf.convertIn(11);

        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 5);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == 8);
        EXPECT_TRUE(mf.getCanonicalForm(mf.add(x,y)) ==
                                          mf.getCanonicalForm(mf.convertIn(4)));
        EXPECT_TRUE(mf.getUnityValue() == mf.getCanonicalForm(mf.convertIn(1)));
        EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalForm(mf.convertIn(0)));
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                  mf.getCanonicalForm(mf.convertIn(static_cast<T>(modulus-1))));
        EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.multiply(y,x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.square(y)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 0)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 1)) == 11);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 2)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 5)) == 7);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 7)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 8)) == 9);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 11)) == 6);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 12)) == 1);
    }

    // try a few tests with the smallest possible modulus
    {
        T modulus = 3;
        M mf(modulus);
        V x = mf.convertIn(1);
        V y = mf.convertIn(2);

        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == 2);
        EXPECT_TRUE(mf.getCanonicalForm(mf.subtract(x,y)) ==
                                          mf.getCanonicalForm(mf.convertIn(2)));
        EXPECT_TRUE(mf.getUnityValue() == mf.getCanonicalForm(mf.convertIn(1)));
        EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalForm(mf.convertIn(0)));
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                  mf.getCanonicalForm(mf.convertIn(static_cast<T>(modulus-1))));
        EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.multiply(y,x)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.square(y)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 0)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 1)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 2)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 3)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 6)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 17)) == 2);
    }

    // try the largest possible modulus
    {
        T modulus = M::max_modulus();
        M mf(modulus);
        V x = mf.convertIn(static_cast<T>(modulus - 1));
        V y = mf.convertIn(2);
    
        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 3);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == modulus - 3);
        EXPECT_TRUE(mf.getCanonicalForm(mf.add(x,y)) ==
                                        mf.getCanonicalForm(mf.convertIn(1)));
        EXPECT_TRUE(mf.getUnityValue() == mf.getCanonicalForm(mf.convertIn(1)));
        EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalForm(mf.convertIn(0)));
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                  mf.getCanonicalForm(mf.convertIn(static_cast<T>(modulus-1))));
        EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == modulus - 2);
        EXPECT_TRUE(mf.convertOut(mf.square(x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 1)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 2)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 10)) == (1024 % modulus));
    }

    // perform a bunch of general checks

    {
        M mf(11);
        T a=5; T b=6;
        test_mf_general_checks(mf, a, b);
        a=0; b=7;
        test_mf_general_checks(mf, a, b);
        a=0; b=0;
        test_mf_general_checks(mf, a, b);
        a=3; b=8;
        test_mf_general_checks(mf, a, b);
        a=2; b=10;
        test_mf_general_checks(mf, a, b);
        a=7; b=9;
        test_mf_general_checks(mf, a, b);
    }
    {
        M mf(M::max_modulus()-2);
        T a=5; T b=6;
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()-1); b=7;
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()/2); b=static_cast<T>(a+3);
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()/2-1); b=static_cast<T>(a+2);
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()-1); b=0;
        test_mf_general_checks(mf, a, b);
        a=2; b=static_cast<T>(mf.getModulus()-2);
        test_mf_general_checks(mf, a, b);
    }
    {
        M mf((M::max_modulus()/4)*2 + 1);
        T a=5; T b=6;
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()-1); b=3;
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()/2); b=static_cast<T>(a+3);
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()/2-1); b=static_cast<T>(a+2);
        test_mf_general_checks(mf, a, b);
        a=static_cast<T>(mf.getModulus()-1); b=0;
        test_mf_general_checks(mf, a, b);
        a=2; b=static_cast<T>(mf.getModulus()-2);
        test_mf_general_checks(mf, a, b);
    }
}


template <template<class> class M>
void test_custom_monty()
{
    namespace mont = hurchalla::montgomery_arithmetic;
    test_MontgomeryForm<mont::MontgomeryForm<uint8_t, M<uint8_t>>>();
    test_MontgomeryForm<mont::MontgomeryForm<uint16_t, M<uint16_t>>>();
    test_MontgomeryForm<mont::MontgomeryForm<uint32_t, M<uint32_t>>>();
    test_MontgomeryForm<mont::MontgomeryForm<uint64_t, M<uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<mont::MontgomeryForm<__uint128_t, M<__uint128_t>>>();
#endif
}


namespace {
    TEST(MontgomeryArithmetic, MontyDefault) {
        namespace ma = hurchalla::modular_arithmetic;
        namespace mont = hurchalla::montgomery_arithmetic;
        test_MontgomeryForm<mont::MontgomeryForm<uint8_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<uint16_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<uint32_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<uint64_t>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_MontgomeryForm<mont::MontgomeryForm<__uint128_t>>();
#endif
        test_MontgomeryForm<mont::MontgomeryForm<int8_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<int16_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<int32_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<int64_t>>();
// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_MontgomeryForm<mont::MontgomeryForm<__int128_t>>();
#endif
    }

    TEST(MontgomeryArithmetic, MontyWrappedStandardMath) {
        namespace mont = hurchalla::montgomery_arithmetic;
        test_custom_monty<mont::MontyWrappedStandardMath>();
    }

    TEST(MontgomeryArithmetic, MontyFullRange) {
        namespace mont = hurchalla::montgomery_arithmetic;
        test_custom_monty<mont::MontyFullRange>();
    }

    TEST(MontgomeryArithmetic, MontyHalfRange) {
        namespace mont = hurchalla::montgomery_arithmetic;
        test_custom_monty<mont::MontyHalfRange>();
    }

    TEST(MontgomeryArithmetic, MontyQuarterRange) {
        namespace mont = hurchalla::montgomery_arithmetic;
        test_custom_monty<mont::MontyQuarterRange>();
    }

    TEST(MontgomeryArithmetic, MontySqrtRange) {
        namespace mont = hurchalla::montgomery_arithmetic;
        test_custom_monty<mont::MontySqrtRange>();
    }
}
