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
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>


template <typename M>
void test_multiply_variants(const M& mf, typename M::MontgomeryValue x,
              typename M::MontgomeryValue y, typename M::T_type expected_result)
{
    namespace hc = hurchalla;
    EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
             mf.template multiply<hc::LowlatencyTag>(x,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
             mf.template multiply<hc::LowuopsTag>(x,y)) == expected_result);

    typename M::MontgomeryValue result;
    bool isZero;
    result = mf.multiplyIsZero(x,y,isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template multiplyIsZero<hc::LowlatencyTag>(x,y,isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template multiplyIsZero<hc::LowuopsTag>(x,y,isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
}

template <typename M>
void test_fmadd_variants(const M& mf, typename M::MontgomeryValue x,
                   typename M::MontgomeryValue y, typename M::CanonicalValue zc,
                   typename M::T_type expected_result)
{
    namespace hc = hurchalla;
    EXPECT_TRUE(mf.convertOut(mf.fmadd(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmadd<hc::LowlatencyTag>(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmadd<hc::LowuopsTag>(x,y,zc)) == expected_result);
}

template <typename M>
void test_fmsub_variants(const M& mf, typename M::MontgomeryValue x,
                   typename M::MontgomeryValue y, typename M::CanonicalValue zc,
                   typename M::T_type expected_result)
{
    namespace hc = hurchalla;
    EXPECT_TRUE(mf.convertOut(mf.fmsub(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmsub<hc::LowlatencyTag>(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmsub<hc::LowuopsTag>(x,y,zc)) == expected_result);
}

template <typename M>
void test_famul_variants(const M& mf, typename M::MontgomeryValue x,
                   typename M::CanonicalValue yc, typename M::MontgomeryValue z,
                   typename M::T_type expected_result)
{
    namespace hc = hurchalla;
    typename M::MontgomeryValue result;
    bool isZero;

    EXPECT_TRUE(mf.convertOut(mf.template
          famul<false>(x,yc,z)) == expected_result);
    EXPECT_TRUE(mf.convertOut(mf.template
          famul<false, hc::LowlatencyTag>(x,yc,z)) == expected_result);
    EXPECT_TRUE(mf.convertOut(mf.template
          famul<false, hc::LowuopsTag>(x,yc,z)) == expected_result);

    result = mf.template famulIsZero<false>(x,yc,z,isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template famulIsZero<false, hc::LowlatencyTag>(x,yc,z,isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template famulIsZero<false, hc::LowuopsTag>(x,yc,z,isZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));

    // This clause is a bit of a hack, since it requires white-box knowledge
    // that modulus < max_modulus()/2 will satisfy the preconditions of all
    // known implementing MontyType famul functions.  New MontyTypes might be
    // added or changed, or some MontyType outside my awareness might be used,
    // all of which could break this assumption.
    if (mf.getModulus() < mf.max_modulus()/2) {
        EXPECT_TRUE(mf.convertOut(mf.template
              famul<true>(x,yc,z)) == expected_result);
        EXPECT_TRUE(mf.convertOut(mf.template
              famul<true, hc::LowlatencyTag>(x,yc,z)) == expected_result);
        EXPECT_TRUE(mf.convertOut(mf.template
              famul<true, hc::LowuopsTag>(x,yc,z)) == expected_result);

        result = mf.template famulIsZero<true>(x,yc,z,isZero);
        EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
        result = mf.template famulIsZero<true,hc::LowlatencyTag>(x,yc,z,isZero);
        EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
        result = mf.template famulIsZero<true, hc::LowuopsTag>(x,yc,z,isZero);
        EXPECT_TRUE(mf.convertOut(result) == expected_result &&
                 isZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    }
}



template <typename M>
void test_mf_general_checks(M& mf, typename M::T_type a, typename M::T_type b,
                                                           typename M::T_type c)
{
    namespace hc = hurchalla;

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

    T reference_sum = hc::modular_addition_prereduced_inputs(a,b,modulus);
    EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == reference_sum);
    EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,y)) ==
                             mf.getCanonicalValue(mf.convertIn(reference_sum)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,yc)) ==
                             mf.getCanonicalValue(mf.convertIn(reference_sum)));

    T diff1 = hc::modular_subtraction_prereduced_inputs(b,a,modulus);
    EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == diff1);
    EXPECT_TRUE(mf.convertOut(mf.subtract(y,xc)) == diff1);
    T diff2 = hc::modular_subtraction_prereduced_inputs(a,b,modulus);
    EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == diff2);
    EXPECT_TRUE(mf.convertOut(mf.subtract(x,yc)) == diff2);
    T us = mf.convertOut(mf.unorderedSubtract(x,y));
    EXPECT_TRUE(us == diff1 || us == diff2);
    us = mf.convertOut(mf.unorderedSubtract(y,x));
    EXPECT_TRUE(us == diff1 || us == diff2);

    EXPECT_TRUE(mf.getUnityValue() == mf.getCanonicalValue(mf.convertIn(1)));
    EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalValue(mf.convertIn(0)));
    EXPECT_TRUE(mf.getNegativeOneValue() ==
              mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));

    T ref_product = hc::modular_multiplication_prereduced_inputs(a,b,modulus);
    test_multiply_variants(mf, x, y, ref_product);
    test_multiply_variants(mf, y, x, ref_product);
    test_fmadd_variants(mf, x, y, zc,
                 hc::modular_addition_prereduced_inputs(ref_product,c,modulus));
    test_fmsub_variants(mf, x, y, zc,
              hc::modular_subtraction_prereduced_inputs(ref_product,c,modulus));
    test_famul_variants(mf, x, yc, z,
          hc::modular_multiplication_prereduced_inputs(
              hc::modular_addition_prereduced_inputs(a,b,modulus), c, modulus));

    T a_squared = hc::modular_multiplication_prereduced_inputs(a,a,modulus);
    test_multiply_variants(mf, x, x, a_squared);
    test_fmadd_variants(mf, x, x, zc,
                   hc::modular_addition_prereduced_inputs(a_squared,c,modulus));
    test_fmsub_variants(mf, x, x, zc,
                hc::modular_subtraction_prereduced_inputs(a_squared,c,modulus));
    test_famul_variants(mf, x, xc, z,
          hc::modular_multiplication_prereduced_inputs(
              hc::modular_addition_prereduced_inputs(a,a,modulus), c, modulus));

    T b_squared = hc::modular_multiplication_prereduced_inputs(b,b,modulus);
    test_multiply_variants(mf, y, y, b_squared);
    test_fmadd_variants(mf, y, y, zc,
                   hc::modular_addition_prereduced_inputs(b_squared,c,modulus));
    test_fmsub_variants(mf, y, y, zc,
                hc::modular_subtraction_prereduced_inputs(b_squared,c,modulus));
    test_famul_variants(mf, y, yc, z,
          hc::modular_multiplication_prereduced_inputs(
              hc::modular_addition_prereduced_inputs(b,b,modulus), c, modulus));

    EXPECT_TRUE(mf.convertOut(mf.pow(y,0)) == 1);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,1)) == b);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,4)) == hc::modular_pow<T>(b,4,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,13)) ==
                                              hc::modular_pow<T>(b,13,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,17)) ==
                                              hc::modular_pow<T>(b,17,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,127)) ==
                                             hc::modular_pow<T>(b,127,modulus));
}


// an example functor to use while testing MontgomeryForm's gcd()
template <typename T>
struct GcdFunctor {
    T operator()(T a, T b)
    {
        static_assert(hurchalla::ut_numeric_limits<T>::is_integer, "");
        static_assert(!hurchalla::ut_numeric_limits<T>::is_signed, "");
        HPBC_PRECONDITION2(a > 0 || b > 0);
        while (a != 0) {
            T tmp = a;
            a = static_cast<T>(b % a);
            b = tmp;
        }
        HPBC_POSTCONDITION2(b > 0);
        return b;
    }
};



template <typename M>
void test_MontgomeryForm()
{
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
        T us = mf.convertOut(mf.unorderedSubtract(x,y));
        EXPECT_TRUE(us == 8 || us == 5);
        us = mf.convertOut(mf.unorderedSubtract(y,x));
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
        EXPECT_TRUE(mf.getCanonicalValue(mf.template
            famulIsZero<false>(one,zero,one,isZero)) == one && isZero == false);
        EXPECT_TRUE(mf.getCanonicalValue(mf.template
            famulIsZero<false>(one,one,zero,isZero)) == zero && isZero == true);
        // This is a bit of a hack, since it requires white-box knowledge that
        // modulus < max_modulus()/2 will satisfy the preconditions of all
        // known implementing MontyType famul functions.  New MontyTypes might
        // be added or changed, or some MontyType outside my awareness might be
        // used, all of which could break this assumption.
        if (mf.getModulus() < mf.max_modulus()/2) {
            EXPECT_TRUE(mf.getCanonicalValue(mf.template
             famulIsZero<true>(one,zero,one,isZero)) == one && isZero == false);
            EXPECT_TRUE(mf.getCanonicalValue(mf.template
             famulIsZero<true>(one,one,zero,isZero)) == zero && isZero == true);
        }

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
        T us = mf.convertOut(mf.unorderedSubtract(x,y));
        EXPECT_TRUE(us == 1 || us == 2);
        us = mf.convertOut(mf.unorderedSubtract(y,x));
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
        T us = mf.convertOut(mf.unorderedSubtract(x,y));
        EXPECT_TRUE(us == 3 || us == modulus - 3);
        us = mf.convertOut(mf.unorderedSubtract(y,x));
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

    // test gcd
    {
        M mf(static_cast<T>(35));
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(28)) == 7);
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(29)) == 1);
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(70)) == 35);
    }
    if (117 <= M::max_modulus())
    {
        M mf(static_cast<T>(117));
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(78)) == 39);
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(26)) == 13);
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(27)) == 9);
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(28)) == 1);
    }
    {
        M mf(static_cast<T>(3));  // smallest possible modulus
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(2)) == 1);
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(0)) == 3);
    }
    {
        T modulus = M::max_modulus();
        ASSERT_TRUE(modulus > 9);
        while (modulus % 3 != 0 || modulus % 2 == 0)
            --modulus;
        M mf(modulus);
        EXPECT_TRUE(
              mf.template gcd_with_modulus<GcdFunctor>(mf.convertIn(12)) == 3);
    }
}


template <template<class> class M>
void test_custom_monty()
{
    namespace hc = hurchalla;
    test_MontgomeryForm<hc::MontgomeryForm<std::uint8_t, M<std::uint8_t>>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint16_t, M<std::uint16_t>>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint32_t, M<std::uint32_t>>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint64_t, M<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<hc::MontgomeryForm<__uint128_t, M<__uint128_t>>>();
#endif
}


namespace {
    TEST(MontgomeryArithmetic, MontyDefault) {
        namespace hc = hurchalla;
        test_MontgomeryForm<hc::MontgomeryForm<std::uint8_t>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::uint16_t>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::uint32_t>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::uint64_t>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_MontgomeryForm<hc::MontgomeryForm<__uint128_t>>();
#endif
        test_MontgomeryForm<hc::MontgomeryForm<std::int8_t>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::int16_t>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::int32_t>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::int64_t>>();
// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_MontgomeryForm<hc::MontgomeryForm<__int128_t>>();
#endif
    }

    TEST(MontgomeryArithmetic, MontyWrappedStandardMath) {
        namespace hc = hurchalla;
        test_custom_monty<hc::detail::MontyWrappedStandardMath>();
    }

    TEST(MontgomeryArithmetic, MontyFullRange) {
        namespace hc = hurchalla;
        test_custom_monty<hc::detail::MontyFullRange>();
    }

    TEST(MontgomeryArithmetic, MontyQuarterRange) {
        namespace hc = hurchalla;
        test_custom_monty<hc::detail::MontyQuarterRange>();
    }

    TEST(MontgomeryArithmetic, MontySqrtRange) {
        namespace hc = hurchalla;
        test_MontgomeryForm<hc::MontgomeryForm<std::uint8_t,
                                  hc::detail::MontySqrtRange<std::uint16_t>>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::uint16_t,
                                  hc::detail::MontySqrtRange<std::uint32_t>>>();
        test_MontgomeryForm<hc::MontgomeryForm<std::uint32_t,
                                  hc::detail::MontySqrtRange<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_MontgomeryForm<hc::MontgomeryForm<std::uint64_t,
                                    hc::detail::MontySqrtRange<__uint128_t>>>();
#endif
    }
}
