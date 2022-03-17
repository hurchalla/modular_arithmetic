// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// Strictly for testing purposes, we'll define HURCHALLA_ALLOW_INLINE_ASM_ALL,
// which will cause MontgomeryForm to use all helper inline asm functions that
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
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontyFullRangeMasked.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace {


namespace hc = ::hurchalla;

// These adapter functions exist because the functions in the modular_arithmetic
// library require unsigned integers, and the tests in this file sometimes need
// signed integers.
namespace testmf_adapters {
    template <typename T>
    T modadd(T a, T b, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_PRECONDITION(0 <= a && a < modulus);
        HPBC_PRECONDITION(0 <= b && b < modulus);
        HPBC_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_addition_prereduced_inputs(
                static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus)));

        HPBC_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }

    template <typename T>
    T modsub(T a, T b, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_PRECONDITION(0 <= a && a < modulus);
        HPBC_PRECONDITION(0 <= b && b < modulus);
        HPBC_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_subtraction_prereduced_inputs(
                static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus)));

        HPBC_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }

    template <typename T>
    T modmul(T a, T b, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_PRECONDITION(0 <= a && a < modulus);
        HPBC_PRECONDITION(0 <= b && b < modulus);
        HPBC_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_multiplication_prereduced_inputs(
                static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus)));

        HPBC_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }

    template <typename T>
    T modpow(T base, T exponent, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_PRECONDITION(base >= 0);
        HPBC_PRECONDITION(exponent >= 0);
        HPBC_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_pow(static_cast<U>(base),
                            static_cast<U>(exponent), static_cast<U>(modulus)));

        HPBC_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }
}



template <typename M>
void test_multiply_variants(const M& mf, typename M::MontgomeryValue x,
         typename M::MontgomeryValue y, typename M::IntegerType expected_result)
{
    EXPECT_TRUE(mf.convertOut(mf.multiply(x,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
             mf.template multiply<hc::LowlatencyTag>(x,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
             mf.template multiply<hc::LowuopsTag>(x,y)) == expected_result);

    typename M::MontgomeryValue result;
    bool resultIsZero;
    result = mf.multiply(x,y,resultIsZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
           resultIsZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template multiply<hc::LowlatencyTag>(x,y,resultIsZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
           resultIsZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
    result = mf.template multiply<hc::LowuopsTag>(x,y,resultIsZero);
    EXPECT_TRUE(mf.convertOut(result) == expected_result &&
           resultIsZero == (mf.getCanonicalValue(result) == mf.getZeroValue()));
}

template <typename M>
void test_fmadd_variants(const M& mf, typename M::MontgomeryValue x,
            typename M::MontgomeryValue y, typename M::CanonicalValue zc,
            typename M::FusingValue zf, typename M::IntegerType expected_result)
{
    EXPECT_TRUE(mf.convertOut(mf.fmadd(x,y,zf)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmadd<hc::LowlatencyTag>(x,y,zf)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmadd<hc::LowuopsTag>(x,y,zf)) == expected_result);

    EXPECT_TRUE(mf.convertOut(mf.fmadd(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmadd<hc::LowlatencyTag>(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmadd<hc::LowuopsTag>(x,y,zc)) == expected_result);
}

template <typename M>
void test_fmsub_variants(const M& mf, typename M::MontgomeryValue x,
            typename M::MontgomeryValue y, typename M::CanonicalValue zc,
            typename M::FusingValue zf, typename M::IntegerType expected_result)
{
    EXPECT_TRUE(mf.convertOut(mf.fmsub(x,y,zf)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmsub<hc::LowlatencyTag>(x,y,zf)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmsub<hc::LowuopsTag>(x,y,zf)) == expected_result);

    EXPECT_TRUE(mf.convertOut(mf.fmsub(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmsub<hc::LowlatencyTag>(x,y,zc)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
          mf.template fmsub<hc::LowuopsTag>(x,y,zc)) == expected_result);
}


template <typename M>
void test_square_variants(const M& mf, typename M::MontgomeryValue x,
                                                  typename M::CanonicalValue zc)
{
    typename M::CanonicalValue answer;

    answer = mf.getCanonicalValue(mf.multiply(x,x));
    EXPECT_TRUE(mf.getCanonicalValue(mf.square(x)) == answer);
    EXPECT_TRUE(mf.getCanonicalValue(mf.template square<hc::LowlatencyTag>(x))
                == answer);
    EXPECT_TRUE(mf.getCanonicalValue(mf.template square<hc::LowuopsTag>(x))
                == answer);

    answer = mf.getCanonicalValue(mf.subtract(mf.multiply(x,x),zc));
    EXPECT_TRUE(mf.getCanonicalValue(mf.fusedSquareSub(x,zc)) == answer);
    EXPECT_TRUE(mf.getCanonicalValue(mf.template
                            fusedSquareSub<hc::LowlatencyTag>(x,zc)) == answer);
    EXPECT_TRUE(mf.getCanonicalValue(mf.template
                            fusedSquareSub<hc::LowuopsTag>(x,zc)) == answer);

    answer = mf.getCanonicalValue(mf.add(mf.multiply(x,x),zc));
    EXPECT_TRUE(mf.getCanonicalValue(mf.fusedSquareAdd(x,zc)) == answer);
    EXPECT_TRUE(mf.getCanonicalValue(mf.template
                            fusedSquareAdd<hc::LowlatencyTag>(x,zc)) == answer);
    EXPECT_TRUE(mf.getCanonicalValue(mf.template
                            fusedSquareAdd<hc::LowuopsTag>(x,zc)) == answer);
}


template <typename M>
void test_mf_general_checks(M& mf, typename M::IntegerType a,
                           typename M::IntegerType b, typename M::IntegerType c)
{
    namespace tma = ::testmf_adapters;

    using T = typename M::IntegerType;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;
    using FV = typename M::FusingValue;
    T modulus = mf.getModulus();
    V x = mf.convertIn(a);
    V y = mf.convertIn(b);
    V z = mf.convertIn(c);
    C xc = mf.getCanonicalValue(x);
    C yc = mf.getCanonicalValue(y);
    C zc = mf.getCanonicalValue(z);
    FV zf = mf.getFusingValue(z);

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

    T reference_sum = tma::modadd(a, b, modulus);
    EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == reference_sum);
    EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == reference_sum);
    EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,y)) ==
                             mf.getCanonicalValue(mf.convertIn(reference_sum)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,yc)) ==
                             mf.getCanonicalValue(mf.convertIn(reference_sum)));

    T diff1 = tma::modsub(b, a, modulus);
    EXPECT_TRUE(mf.convertOut(mf.subtract(y,x)) == diff1);
    EXPECT_TRUE(mf.convertOut(mf.subtract(y,xc)) == diff1);
    T diff2 = tma::modsub(a, b, modulus);
    EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == diff2);
    EXPECT_TRUE(mf.convertOut(mf.subtract(x,yc)) == diff2);
    T us = mf.convertOut(mf.unorderedSubtract(x,y));
    EXPECT_TRUE(us == diff1 || us == diff2);
    us = mf.convertOut(mf.unorderedSubtract(y,x));
    EXPECT_TRUE(us == diff1 || us == diff2);

    EXPECT_TRUE(mf.getUnityValue() == mf.getCanonicalValue(mf.convertIn(1)));
    EXPECT_TRUE(mf.getZeroValue() == mf.getCanonicalValue(mf.convertIn(0)));
    EXPECT_TRUE(modulus > 0);
    EXPECT_TRUE(mf.getNegativeOneValue() ==
              mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));

    T ref_product = tma::modmul(a, b, modulus);
    test_multiply_variants(mf, x, y, ref_product);
    test_multiply_variants(mf, y, x, ref_product);
    test_fmadd_variants(mf, x, y, zc, zf, tma::modadd(ref_product,c,modulus));
    test_fmsub_variants(mf, x, y, zc, zf, tma::modsub(ref_product,c,modulus));

    T a_squared = tma::modmul(a, a, modulus);
    test_multiply_variants(mf, x, x, a_squared);
    test_fmadd_variants(mf, x, x, zc, zf, tma::modadd(a_squared,c,modulus));
    test_fmsub_variants(mf, x, x, zc, zf, tma::modsub(a_squared,c,modulus));

    T b_squared = tma::modmul(b, b, modulus);
    test_multiply_variants(mf, y, y, b_squared);
    test_fmadd_variants(mf, y, y, zc, zf, tma::modadd(b_squared,c,modulus));
    test_fmsub_variants(mf, y, y, zc, zf, tma::modsub(b_squared,c,modulus));

    test_square_variants(mf, x, zc);
    test_square_variants(mf, y, zc);

    EXPECT_TRUE(mf.convertOut(mf.pow(y,0)) == 1);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,1)) == b);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,4)) == tma::modpow<T>(b,4,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,13)) == tma::modpow<T>(b,13,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,17)) == tma::modpow<T>(b,17,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,127)) == tma::modpow<T>(b,127,modulus));
#ifdef __GNUC__
#  pragma GCC diagnostic push
// old versions of gcc and clang give unnecessary warnings about single braced
// initialization lists with std::array (newer versions fixed this).
#  pragma GCC diagnostic ignored "-Wmissing-braces"
#endif
    std::array<V,3> mv_base = { x, y, z };
#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif
    // Do just a simple test of pow()'s array form template function -
    // it's tested more thoroughly in test_montgomery_pow.cpp.
    std::array<V,3> mv_res = mf.pow(mv_base, 19);
    EXPECT_TRUE(mf.convertOut(mv_res[0]) == tma::modpow<T>(a,19,modulus));
    EXPECT_TRUE(mf.convertOut(mv_res[1]) == tma::modpow<T>(b,19,modulus));
    EXPECT_TRUE(mf.convertOut(mv_res[2]) == tma::modpow<T>(c,19,modulus));
}



// an example functor to use while testing MontgomeryForm's gcd()
struct GcdFunctor {
    template <typename T>
    T operator()(T a, T b) const
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        static_assert(!hc::ut_numeric_limits<T>::is_signed, "");
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
    using T = typename M::IntegerType;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;
    using FV = typename M::FusingValue;

    // Try a basic test case first that is valid for all possible Monty types
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
        EXPECT_TRUE(modulus > 0);
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                 mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));
        test_multiply_variants(mf, x, y, 1);
        test_multiply_variants(mf, y, x, 1);
        test_multiply_variants(mf, y, y, 4);

        V z = mf.convertIn(9);
        C zc = mf.getCanonicalValue(z);
        FV zf = mf.getFusingValue(z);
        test_fmadd_variants(mf, x, y, zc, zf, 10);
        test_fmadd_variants(mf, x, x, zc, zf, 6);
        test_fmadd_variants(mf, y, y, zc, zf, 0);
        test_fmsub_variants(mf, x, y, zc, zf, 5);
        test_fmsub_variants(mf, x, x, zc, zf, 1);
        test_fmsub_variants(mf, y, y, zc, zf, 8);

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

        // test whether the multiply overload correctly sets resultIsZero 
        C zero = mf.getZeroValue();
        C one = mf.getUnityValue();
        bool resultIsZero;
        EXPECT_TRUE(mf.getCanonicalValue(mf.multiply(one,one,resultIsZero)) ==
                    one && resultIsZero == false);
        EXPECT_TRUE(mf.getCanonicalValue(mf.multiply(one,zero,resultIsZero)) ==
                    zero && resultIsZero == true);
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
        EXPECT_TRUE(modulus > 0);
        EXPECT_TRUE(mf.getNegativeOneValue() ==
                 mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-1))));
        test_multiply_variants(mf, x, y, 2);
        test_multiply_variants(mf, y, x, 2);
        test_multiply_variants(mf, y, y, 1);

        V z = mf.convertIn(1);
        C zc = mf.getCanonicalValue(z);
        FV zf = mf.getFusingValue(z);
        test_fmadd_variants(mf, x, y, zc, zf, 0);
        test_fmadd_variants(mf, x, x, zc, zf, 2);
        test_fmadd_variants(mf, y, y, zc, zf, 2);
        test_fmsub_variants(mf, x, y, zc, zf, 1);
        test_fmsub_variants(mf, x, x, zc, zf, 0);
        test_fmsub_variants(mf, y, y, zc, zf, 0);

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
        EXPECT_TRUE(modulus > 4);

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
        FV zf = mf.getFusingValue(z);
        test_fmadd_variants(mf, x, y, zc, zf, static_cast<T>(modulus - 1));
        test_fmadd_variants(mf, x, x, zc, zf, 2);
        test_fmsub_variants(mf, x, y, zc, zf, static_cast<T>(modulus - 3));
        test_fmsub_variants(mf, x, x, zc, zf, 0);

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
        EXPECT_TRUE(M::max_modulus() >= 5);
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
        EXPECT_TRUE(M::max_modulus() >= 5);
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
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(28), GcdFunctor()) == 7);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(29), GcdFunctor()) == 1);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(70), GcdFunctor()) == 35);
    }
    if (117 <= M::max_modulus())
    {
        M mf(static_cast<T>(117));
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(78), GcdFunctor()) == 39);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(26), GcdFunctor()) == 13);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(27), GcdFunctor()) == 9);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(28), GcdFunctor()) == 1);
    }
    {
        M mf(static_cast<T>(3));  // smallest possible modulus
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(2), GcdFunctor()) == 1);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(0), GcdFunctor()) == 3);
    }
    {
        T modulus = M::max_modulus();
        ASSERT_TRUE(modulus > 9);
        while (modulus % 3 != 0 || modulus % 2 == 0)
            --modulus;
        M mf(modulus);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(12), GcdFunctor()) == 3);
    }
}


template <template<class> class M>
void test_custom_monty()
{
    test_MontgomeryForm<hc::MontgomeryForm<std::uint8_t, M<std::uint8_t>>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint16_t, M<std::uint16_t>>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint32_t, M<std::uint32_t>>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint64_t, M<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<hc::MontgomeryForm<__uint128_t, M<__uint128_t>>>();
#endif
}



TEST(MontgomeryArithmetic, MontyDefault) {
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

// check that MontgomeryDefault uses MontyHalfRange when appropriate
#if HURCHALLA_TARGET_BIT_WIDTH == 32
    static_assert(std::is_same<
                        hc::detail::MontgomeryDefault<std::int32_t>::type,
                        hc::detail::MontyHalfRange<std::uint32_t>
                  >::value, "");
#elif HURCHALLA_TARGET_BIT_WIDTH == 64
    static_assert(std::is_same<
                        hc::detail::MontgomeryDefault<std::int64_t>::type,
                        hc::detail::MontyHalfRange<std::uint64_t>
                  >::value, "");
#endif
}

TEST(MontgomeryArithmetic, MontyWrappedStandardMath) {
    test_custom_monty<hc::detail::MontyWrappedStandardMath>();
}

TEST(MontgomeryArithmetic, MontyFullRange) {
    test_custom_monty<hc::detail::MontyFullRange>();
}

TEST(MontgomeryArithmetic, MontyHalfRange) {
    test_custom_monty<hc::detail::MontyHalfRange>();
}

TEST(MontgomeryArithmetic, MontyFullRangeMasked) {
    test_custom_monty<hc::detail::MontyFullRangeMasked>();
}

TEST(MontgomeryArithmetic, MontyQuarterRange) {
    test_custom_monty<hc::detail::MontyQuarterRange>();
}


} // end unnamed namespace
