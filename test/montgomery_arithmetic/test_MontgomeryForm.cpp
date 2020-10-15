// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


// For more complete testing we'll define HURCHALLA_ALLOW_INLINE_ASM_ALL,
// which will cause MontgomeryForm to use any helper inline asm functions that
// are available.  Internally, these inline asm functions will also call their
// corresponding generic template helper functions inside a postcondition, in
// order to make sure that the asm result is correct.  Of course postcondition
// checks must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to undefine NDEBUG, which is why we undef
// NDEBUG here too.
#undef HURCHALLA_ALLOW_INLINE_ASM_ALL
#define HURCHALLA_ALLOW_INLINE_ASM_ALL 1
#undef NDEBUG


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySixthRange.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "gtest/gtest.h"
#include <cstdint>


template <typename M>
void test_multiply_variants(const M& mf, typename M::MontgomeryValue x,
              typename M::MontgomeryValue y, typename M::T_type expected_result)
{
    namespace ma = hurchalla::montgomery_arithmetic;
    EXPECT_TRUE(mf.convertOut(mf.multiply(x, y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
      mf.template multiply<ma::LowlatencyTag>(x,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
      mf.template multiply<ma::LowuopsTag>(x,y)) == expected_result);

    typename M::MontgomeryValue result;
    bool isZero;
    result = mf.multiplyIsZero(x, y, isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template multiplyIsZero<ma::LowlatencyTag>(x, y, isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template multiplyIsZero<ma::LowuopsTag>(x, y, isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
}

template <typename M>
void test_fmadd_variants(const M& mf, typename M::MontgomeryValue x,
                   typename M::MontgomeryValue y, typename M::CanonicalValue zc,
                   typename M::T_type expected_result)
{
    namespace ma = hurchalla::montgomery_arithmetic;
    EXPECT_TRUE(mf.convertOut(mf.fmadd(x, y, zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
       mf.template fmadd<ma::LowlatencyTag>(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
      mf.template fmadd<ma::LowuopsTag>(x,y,zc))==expected_result);
}

template <typename M>
void test_fmsub_variants(const M& mf, typename M::MontgomeryValue x,
                   typename M::MontgomeryValue y, typename M::CanonicalValue zc,
                   typename M::T_type expected_result)
{
    namespace ma = hurchalla::montgomery_arithmetic;
    EXPECT_TRUE(mf.convertOut(mf.fmsub(x, y, zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
       mf.template fmsub<ma::LowlatencyTag>(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
      mf.template fmsub<ma::LowuopsTag>(x,y,zc))==expected_result);
}

template <typename M>
void test_famul_variants(const M& mf, typename M::MontgomeryValue x,
                   typename M::CanonicalValue yc, typename M::MontgomeryValue z,
                   typename M::T_type expected_result)
{
    namespace ma = hurchalla::montgomery_arithmetic;
    EXPECT_TRUE(mf.convertOut(mf.famul(x, yc, z)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
       mf.template famul<ma::LowlatencyTag>(x,yc,z)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
      mf.template famul<ma::LowuopsTag>(x,yc,z))==expected_result);

    typename M::MontgomeryValue result;
    bool isZero;
    result = mf.famulIsZero(x, yc, z, isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template famulIsZero<ma::LowlatencyTag>(x, yc, z, isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template famulIsZero<ma::LowuopsTag>(x, yc, z, isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
}



template <typename M>
void test_mf_general_checks(M& mf, typename M::T_type a, typename M::T_type b,
                                                           typename M::T_type c)
{
    namespace ma = hurchalla::modular_arithmetic;

    using T = typename M::T_type;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;
    T modulus = mf.getModulus();
    V x = mf.convertIn(a);
    V y = mf.convertIn(b);
    V z = mf.convertIn(c);
    C xc = mf.getCanonicalValue(x);
    C yc = mf.getCanonicalValue(y);
    C zc = mf.getCanonicalValue(z);

    EXPECT_TRUE(mf.getCanonicalValue(mf.negate(x)) ==
                mf.getCanonicalValue(mf.subtract(mf.getZeroValue(), x)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.negate(y)) ==
                mf.getCanonicalValue(mf.subtract(mf.getZeroValue(), y)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.negate(z)) ==
                mf.getCanonicalValue(mf.subtract(mf.getZeroValue(), z)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.negate(xc)) ==
                mf.getCanonicalValue(mf.subtract(mf.getZeroValue(), xc)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.negate(yc)) ==
                mf.getCanonicalValue(mf.subtract(mf.getZeroValue(), yc)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.negate(zc)) ==
                mf.getCanonicalValue(mf.subtract(mf.getZeroValue(), zc)));

    T reference_sum = ma::modular_addition_prereduced_inputs(a,b,modulus);
    EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == reference_sum);
    EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,y)) ==
                             mf.getCanonicalValue(mf.convertIn(reference_sum)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,yc)) ==
                             mf.getCanonicalValue(mf.convertIn(reference_sum)));

    T diff1 = ma::modular_subtraction_prereduced_inputs(b,a,modulus);
    EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == diff1);
    EXPECT_TRUE(mf.convertOut(mf.subtract(y,xc)) == diff1);
    T diff2 = ma::modular_subtraction_prereduced_inputs(a,b,modulus);
    EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == diff2);
    EXPECT_TRUE(mf.convertOut(mf.subtract(x,yc)) == diff2);
    T us = mf.convertOut(mf.unordered_subtract(x,y));
    EXPECT_TRUE(us == diff1 || us == diff2);
    us = mf.convertOut(mf.unordered_subtract(y,x));
    EXPECT_TRUE(us == diff1 || us == diff2);

    EXPECT_TRUE(mf.getUnityValue() == mf.getCanonicalValue(mf.convertIn(1)));
    EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalValue(mf.convertIn(0)));
    EXPECT_TRUE(mf.getNegativeOneValue() ==
              mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));

    T ref_product = ma::modular_multiplication_prereduced_inputs(a,b,modulus);
    test_multiply_variants(mf, x, y, ref_product);
    test_multiply_variants(mf, y, x, ref_product);
    test_fmadd_variants(mf, x, y, zc,
                 ma::modular_addition_prereduced_inputs(ref_product,c,modulus));
    test_fmsub_variants(mf, x, y, zc,
              ma::modular_subtraction_prereduced_inputs(ref_product,c,modulus));
    test_famul_variants(mf, x, yc, z,
          ma::modular_multiplication_prereduced_inputs(
            ma::modular_addition_prereduced_inputs(a, b, modulus), c, modulus));

    T a_squared = ma::modular_multiplication_prereduced_inputs(a,a,modulus);
    test_multiply_variants(mf, x, x, a_squared);
    test_fmadd_variants(mf, x, x, zc,
                   ma::modular_addition_prereduced_inputs(a_squared,c,modulus));
    test_fmsub_variants(mf, x, x, zc,
                ma::modular_subtraction_prereduced_inputs(a_squared,c,modulus));
    test_famul_variants(mf, x, xc, z,
          ma::modular_multiplication_prereduced_inputs(
            ma::modular_addition_prereduced_inputs(a, a, modulus), c, modulus));

    T b_squared = ma::modular_multiplication_prereduced_inputs(b,b,modulus);
    test_multiply_variants(mf, y, y, b_squared);
    test_fmadd_variants(mf, y, y, zc,
                   ma::modular_addition_prereduced_inputs(b_squared,c,modulus));
    test_fmsub_variants(mf, y, y, zc,
                ma::modular_subtraction_prereduced_inputs(b_squared,c,modulus));
    test_famul_variants(mf, y, yc, z,
          ma::modular_multiplication_prereduced_inputs(
            ma::modular_addition_prereduced_inputs(b, b, modulus), c, modulus));

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
    namespace ma = hurchalla::montgomery_arithmetic;
    using T = typename M::T_type;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;

    // Try a basic test case first that is valid for all possible Monty types,
    // even M == MontySqrtRange<std::uint8_t>.
    {
        T modulus = 13;
        M mf(modulus);
        V x = mf.convertIn(6);
        V y = mf.convertIn(11);

        C xc = mf.getCanonicalValue(x);
        C yc = mf.getCanonicalValue(y);

        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 5);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == 8);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,xc)) == 5);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,yc)) == 8);
        T us = mf.convertOut(mf.unordered_subtract(x,y));
        EXPECT_TRUE(us == 8 || us == 5);
        us = mf.convertOut(mf.unordered_subtract(y,x));
        EXPECT_TRUE(us == 8 || us == 5);
        EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,y)) ==
                                         mf.getCanonicalValue(mf.convertIn(4)));
        EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,yc)) ==
                                         mf.getCanonicalValue(mf.convertIn(4)));
        EXPECT_TRUE(mf.getUnityValue()== mf.getCanonicalValue(mf.convertIn(1)));
        EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalValue(mf.convertIn(0)));
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                 mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));
        test_multiply_variants(mf, x, y, 1);
        test_multiply_variants(mf, y, x, 1);
        test_multiply_variants(mf, y, y, 4);

        V z = mf.convertIn(9);
        C zc = mf.getCanonicalValue(z);
        test_fmadd_variants(mf, x, y, zc, 10);
        test_fmadd_variants(mf, x, x, zc, 6);
        test_fmadd_variants(mf, y, y, zc, 0);
        test_fmsub_variants(mf, x, y, zc, 5);
        test_fmsub_variants(mf, x, x, zc, 1);
        test_fmsub_variants(mf, y, y, zc, 8);
        test_famul_variants(mf, x, yc, z, 10);
        test_famul_variants(mf, x, xc, z, 4);
        test_famul_variants(mf, y, yc, z, 3);

        EXPECT_TRUE(mf.convertOut(mf.pow(y, 0)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 1)) == 11);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 2)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 5)) == 7);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 7)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 8)) == 9);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 11)) == 6);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 12)) == 1);

        EXPECT_TRUE(mf.convertOut(mf.negate(x)) == 7);
        EXPECT_TRUE(mf.convertOut(mf.negate(xc)) == 7);
        EXPECT_TRUE(mf.convertOut(mf.negate(y)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.negate(yc)) == 2);

        // test to see if famulIsZero and multiplyIsZero set isZero correctly
        C zero = mf.getZeroValue();
        C one = mf.getUnityValue();
        bool isZero;
        EXPECT_TRUE(mf.getCanonicalValue(mf.famulIsZero(one,zero,one,isZero)) ==
                    one && isZero == false);
        EXPECT_TRUE(mf.getCanonicalValue(mf.famulIsZero(one,one,zero,isZero)) ==
                    zero && isZero == true);

        EXPECT_TRUE(mf.getCanonicalValue(mf.multiplyIsZero(one,one,isZero)) ==
                    one && isZero == false);
        EXPECT_TRUE(mf.getCanonicalValue(mf.multiplyIsZero(one,zero,isZero)) ==
                    zero && isZero == true);
    }

    // try tests with the smallest possible modulus
    {
        T modulus = 3;
        M mf(modulus);
        V x = mf.convertIn(1);
        V y = mf.convertIn(2);

        C xc = mf.getCanonicalValue(x);
        C yc = mf.getCanonicalValue(y);

        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,xc)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,yc)) == 2);
        EXPECT_TRUE(mf.getCanonicalValue(mf.subtract(x,y)) ==
                                         mf.getCanonicalValue(mf.convertIn(2)));
        T us = mf.convertOut(mf.unordered_subtract(x,y));
        EXPECT_TRUE(us == 1 || us == 2);
        us = mf.convertOut(mf.unordered_subtract(y,x));
        EXPECT_TRUE(us == 1 || us == 2);
        EXPECT_TRUE(mf.getUnityValue()== mf.getCanonicalValue(mf.convertIn(1)));
        EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalValue(mf.convertIn(0)));
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                 mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));
        test_multiply_variants(mf, x, y, 2);
        test_multiply_variants(mf, y, x, 2);
        test_multiply_variants(mf, y, y, 1);

        V z = mf.convertIn(1);
        C zc = mf.getCanonicalValue(z);
        test_fmadd_variants(mf, x, y, zc, 0);
        test_fmadd_variants(mf, x, x, zc, 2);
        test_fmadd_variants(mf, y, y, zc, 2);
        test_fmsub_variants(mf, x, y, zc, 1);
        test_fmsub_variants(mf, x, x, zc, 0);
        test_fmsub_variants(mf, y, y, zc, 0);
        test_famul_variants(mf, x, yc, z, 0);
        test_famul_variants(mf, x, xc, z, 2);
        test_famul_variants(mf, y, yc, z, 1);

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

        C xc = mf.getCanonicalValue(x);
        C yc = mf.getCanonicalValue(y);
    
        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == 3);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == modulus - 3);
        EXPECT_TRUE(mf.convertOut(mf.subtract(y,xc)) == 3);
        EXPECT_TRUE(mf.convertOut(mf.subtract(x,yc)) == modulus - 3);
        EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,y)) ==
                                        mf.getCanonicalValue(mf.convertIn(1)));
        T us = mf.convertOut(mf.unordered_subtract(x,y));
        EXPECT_TRUE(us == 3 || us == modulus - 3);
        us = mf.convertOut(mf.unordered_subtract(y,x));
        EXPECT_TRUE(us == 3 || us == modulus - 3);
        EXPECT_TRUE(mf.getUnityValue()== mf.getCanonicalValue(mf.convertIn(1)));
        EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalValue(mf.convertIn(0)));
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                 mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));
        test_multiply_variants(mf, x, y, static_cast<T>(modulus - 2));
        test_multiply_variants(mf, x, x, 1);

        V z = mf.convertIn(1);
        C zc = mf.getCanonicalValue(z);
        test_fmadd_variants(mf, x, y, zc, static_cast<T>(modulus - 1));
        test_fmadd_variants(mf, x, x, zc, 2);
        test_fmsub_variants(mf, x, y, zc, static_cast<T>(modulus - 3));
        test_fmsub_variants(mf, x, x, zc, 0);
        test_famul_variants(mf, x, yc, z, 1);
        test_famul_variants(mf, y, xc, z, 1);
        test_famul_variants(mf, x, xc, z, static_cast<T>(modulus - 2));

        EXPECT_TRUE(mf.convertOut(mf.pow(y, 1)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 2)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.pow(y, 10)) == (1024 % modulus));
    }

    // perform a bunch of general checks

    {
        M mf(11);
        T c = 1;
        T a=5; T b=6;
        test_mf_general_checks(mf, a, b, c);
        a=0; b=7;
        test_mf_general_checks(mf, a, b, c);
        a=10; b=0;
        test_mf_general_checks(mf, a, b, c);
        a=0; b=0;
        test_mf_general_checks(mf, a, b, c);
        a=3; b=8;
        test_mf_general_checks(mf, a, b, c);
        a=2; b=10;
        test_mf_general_checks(mf, a, b, c);
        a=7; b=9;
        test_mf_general_checks(mf, a, b, c);
    }
    {
        M mf(M::max_modulus()-2);
        T c = M::max_modulus()-3;
        T a=5; T b=6;
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()-1); b=7;
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()/2); b=static_cast<T>(a+3);
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()/2-1); b=static_cast<T>(a+2);
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()-1); b=0;
        test_mf_general_checks(mf, a, b, c);
        a=2; b=static_cast<T>(mf.getModulus()-2);
        test_mf_general_checks(mf, a, b, c);
    }
    {
        M mf((M::max_modulus()/4)*2 + 1);
        T c = 0;
        T a=5; T b=6;
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()-1); b=3;
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()/2); b=static_cast<T>(a+3);
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()/2-1); b=static_cast<T>(a+2);
        test_mf_general_checks(mf, a, b, c);
        a=static_cast<T>(mf.getModulus()-1); b=0;
        test_mf_general_checks(mf, a, b, c);
        a=2; b=static_cast<T>(mf.getModulus()-2);
        test_mf_general_checks(mf, a, b, c);
    }
}


template <template<class> class M>
void test_custom_monty()
{
    namespace mont = hurchalla::montgomery_arithmetic;
    test_MontgomeryForm<mont::MontgomeryForm<std::uint8_t, M<std::uint8_t>>>();
    test_MontgomeryForm<mont::MontgomeryForm<std::uint16_t,M<std::uint16_t>>>();
    test_MontgomeryForm<mont::MontgomeryForm<std::uint32_t,M<std::uint32_t>>>();
    test_MontgomeryForm<mont::MontgomeryForm<std::uint64_t,M<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<mont::MontgomeryForm<__uint128_t, M<__uint128_t>>>();
#endif
}


namespace {
    TEST(MontgomeryArithmetic, MontyDefault) {
        namespace ma = hurchalla::modular_arithmetic;
        namespace mont = hurchalla::montgomery_arithmetic;
        test_MontgomeryForm<mont::MontgomeryForm<std::uint8_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::uint16_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::uint32_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::uint64_t>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_MontgomeryForm<mont::MontgomeryForm<__uint128_t>>();
#endif
        test_MontgomeryForm<mont::MontgomeryForm<std::int8_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::int16_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::int32_t>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::int64_t>>();
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

    TEST(MontgomeryArithmetic, MontySixthRange) {
        namespace mont = hurchalla::montgomery_arithmetic;
        test_custom_monty<mont::MontySixthRange>();
    }

    TEST(MontgomeryArithmetic, MontySqrtRange) {
        namespace mont = hurchalla::montgomery_arithmetic;
        test_MontgomeryForm<mont::MontgomeryForm<std::uint8_t,
                                        mont::MontySqrtRange<std::uint16_t>>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::uint16_t,
                                        mont::MontySqrtRange<std::uint32_t>>>();
        test_MontgomeryForm<mont::MontgomeryForm<std::uint32_t,
                                        mont::MontySqrtRange<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_MontgomeryForm<mont::MontgomeryForm<std::uint64_t,
                                          mont::MontySqrtRange<__uint128_t>>>();
#endif
    }
}
