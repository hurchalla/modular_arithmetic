// Copyright (c) 2020-2024 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_TEST_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_TEST_MONTGOMERY_FORM_H_INCLUDED


// Strictly for testing purposes, we'll define HURCHALLA_ALLOW_INLINE_ASM_ALL,
// which will cause MontgomeryForm to use all helper inline asm functions that
// are available.  Internally, these inline asm functions will also call their
// corresponding generic template helper functions inside a postcondition, in
// order to make sure that the asm result is correct.  Of course postcondition

// checks must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to define HURCHALLA_CLOCKWORK_ENABLE_ASSERTS,
// which is why we do so here.  This is all strictly for testing purposes.
#undef HURCHALLA_ALLOW_INLINE_ASM_ALL
#define HURCHALLA_ALLOW_INLINE_ASM_ALL 1

#ifndef HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#  define HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#endif


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>
#include <array>
#include <memory>


// We define the following macro in order to be able to specify a variadic
// template argument for ConcreteMontgomeryForm (if it's used) that matches the
// array pow() sizes that will (or that could) be tested by this file.
//
//#define TESTABLE_ARRAY_POW_SIZES() 1, 2, 3, 4, 5
// for now we'll only test the array version of pow with array size 3
#define TESTABLE_ARRAY_POW_SIZES() 3


// These adapter functions exist because the functions in the modular_arithmetic
// library require unsigned integers, and the tests in this file sometimes need
// signed integers.
namespace testmf_adapters {
    namespace hc = ::hurchalla;
    template <typename T>
    T modadd(T a, T b, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_CLOCKWORK_PRECONDITION(0 <= a && a < modulus);
        HPBC_CLOCKWORK_PRECONDITION(0 <= b && b < modulus);
        HPBC_CLOCKWORK_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_addition_prereduced_inputs(
                static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus)));

        HPBC_CLOCKWORK_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }

    template <typename T>
    T modsub(T a, T b, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_CLOCKWORK_PRECONDITION(0 <= a && a < modulus);
        HPBC_CLOCKWORK_PRECONDITION(0 <= b && b < modulus);
        HPBC_CLOCKWORK_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_subtraction_prereduced_inputs(
                static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus)));

        HPBC_CLOCKWORK_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }

    template <typename T>
    T modmul(T a, T b, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_CLOCKWORK_PRECONDITION(0 <= a && a < modulus);
        HPBC_CLOCKWORK_PRECONDITION(0 <= b && b < modulus);
        HPBC_CLOCKWORK_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_multiplication_prereduced_inputs(
                static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus)));

        HPBC_CLOCKWORK_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }

    template <typename T>
    T modpow(T base, T exponent, T modulus)
    {
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        HPBC_CLOCKWORK_PRECONDITION(base >= 0);
        HPBC_CLOCKWORK_PRECONDITION(exponent >= 0);
        HPBC_CLOCKWORK_PRECONDITION(modulus > 1);

        using U = typename hc::extensible_make_unsigned<T>::type;
        T result = static_cast<T>(hc::modular_pow(static_cast<U>(base),
                            static_cast<U>(exponent), static_cast<U>(modulus)));

        HPBC_CLOCKWORK_POSTCONDITION(0 <= result && result < modulus);
        return result;
    }
}


// functor that helps test MontgomeryForm's gcd()
struct GcdFunctor {
    template <typename T>
    T operator()(T a, T b) const
    {
        namespace hc = ::hurchalla;
        static_assert(hc::ut_numeric_limits<T>::is_integer, "");
        static_assert(!hc::ut_numeric_limits<T>::is_signed, "");
        HPBC_CLOCKWORK_PRECONDITION2(a > 0 || b > 0);
        while (a != 0) {
            T tmp = a;
            a = static_cast<T>(b % a);
            b = tmp;
        }
        HPBC_CLOCKWORK_POSTCONDITION2(b > 0);
        return b;
    }
};




namespace {

template <typename M>
void test_subtract_variants(const M& mf, typename M::MontgomeryValue x,
         typename M::MontgomeryValue y, typename M::IntegerType expected_result)
{
    namespace hc = ::hurchalla;
    using C = typename M::CanonicalValue;
    C cx = mf.getCanonicalValue(x);
    C cy = mf.getCanonicalValue(y);

    EXPECT_TRUE(mf.convertOut(mf.subtract(x,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowlatencyTag>(x,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowuopsTag>(x,y)) == expected_result);

    EXPECT_TRUE(mf.convertOut(mf.subtract(cx,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowlatencyTag>(cx,y)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowuopsTag>(cx,y)) == expected_result);

    EXPECT_TRUE(mf.convertOut(mf.subtract(x,cy)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowlatencyTag>(x,cy)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowuopsTag>(x,cy)) == expected_result);

    EXPECT_TRUE(mf.convertOut(mf.subtract(cx,cy)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowlatencyTag>(cx,cy)) == expected_result);
    EXPECT_TRUE(mf.convertOut(
            mf.template subtract<hc::LowuopsTag>(cx,cy)) == expected_result);
}

template <typename M>
void test_multiply_variants(const M& mf, typename M::MontgomeryValue x,
         typename M::MontgomeryValue y, typename M::IntegerType expected_result)
{
    namespace hc = ::hurchalla;
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
    namespace hc = ::hurchalla;
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
    namespace hc = ::hurchalla;
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
    namespace hc = ::hurchalla;
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
void test_remainder(const M& mf)
{
    using T = typename M::IntegerType;
    namespace hc = ::hurchalla;

    T max = hc::ut_numeric_limits<T>::max();
    T mid = static_cast<T>(max/2);
    T modulus = mf.getModulus();

    EXPECT_TRUE(mf.remainder(0) == (0 % modulus));
    EXPECT_TRUE(mf.remainder(1) == (1 % modulus));
    EXPECT_TRUE(mf.remainder(2) == (2 % modulus));
    EXPECT_TRUE(mf.remainder(static_cast<T>(max-0)) == ((max-0) % modulus));
    EXPECT_TRUE(mf.remainder(static_cast<T>(max-1)) == ((max-1) % modulus));
    EXPECT_TRUE(mf.remainder(static_cast<T>(max-2)) == ((max-2) % modulus));
    EXPECT_TRUE(mf.remainder(static_cast<T>(mid-1)) == ((mid-1) % modulus));
    EXPECT_TRUE(mf.remainder(static_cast<T>(mid-0)) == ((mid-0) % modulus));
    EXPECT_TRUE(mf.remainder(static_cast<T>(mid+1)) == ((mid+1) % modulus));
}

template <typename M>
void test_single_inverse(const M& mf, typename M::IntegerType a)
{
    namespace hc = ::hurchalla;
    using T = typename M::IntegerType;
    using U = typename hc::extensible_make_unsigned<T>::type;

    U n = static_cast<U>(mf.getModulus());
    U gcd;  // ignored
    auto answer = hc::modular_multiplicative_inverse(static_cast<U>(a), n, gcd);
    U val = static_cast<U>(mf.convertOut(mf.inverse(mf.convertIn(a))));
    EXPECT_TRUE(val == answer);
}

template <typename M>
void test_inverse(const M& mf)
{
    using T = typename M::IntegerType;
    T max = ::hurchalla::ut_numeric_limits<T>::max();
    T mid = static_cast<T>(max/2);
    T modulus = mf.getModulus();
    test_single_inverse(mf, static_cast<T>(0));
    test_single_inverse(mf, static_cast<T>(1));
    test_single_inverse(mf, static_cast<T>(2));
    test_single_inverse(mf, static_cast<T>(max-0));
    test_single_inverse(mf, static_cast<T>(max-1));
    test_single_inverse(mf, static_cast<T>(mid-0));
    test_single_inverse(mf, static_cast<T>(mid-1));
    test_single_inverse(mf, static_cast<T>(modulus-1));
    test_single_inverse(mf, static_cast<T>(modulus-2));
    test_single_inverse(mf, static_cast<T>(modulus/2));
    test_single_inverse(mf, static_cast<T>((modulus/2) - 1));
}

template <typename M>
void test_divideBySmallPowerOf2_for_dividend(const M& mf,
                                                typename M::IntegerType a)
{
    namespace hc = ::hurchalla;
    using T = typename M::IntegerType;
    using U = typename hc::extensible_make_unsigned<T>::type;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;

    U n = static_cast<U>(mf.getModulus());
    ASSERT_TRUE(n % 2 == 1);
    U gcd;  // ignored
    U inv2 = hc::modular_multiplicative_inverse(static_cast<U>(2), n, gcd);
    ASSERT_TRUE(inv2 != 0);
    V mont_inv2 = mf.convertIn(static_cast<T>(inv2));

    C cx = mf.getCanonicalValue(mf.convertIn(a));
    C answer = cx;
    for (int i=0; i<8; ++i) {
        C val = mf.getCanonicalValue(mf.divideBySmallPowerOf2(cx, i));
        EXPECT_TRUE(val == answer);
        answer = mf.getCanonicalValue(mf.multiply(answer, mont_inv2));
    }
}

template <typename M>
void test_divideBySmallPowerOf2(const M& mf)
{
    using T = typename M::IntegerType;
    T max = ::hurchalla::ut_numeric_limits<T>::max();
    T mid = static_cast<T>(max/2);
    T modulus = mf.getModulus();
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(0));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(1));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(2));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(max-0));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(max-1));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(mid-0));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(mid-1));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(modulus-1));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(modulus-2));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>(modulus/2));
    test_divideBySmallPowerOf2_for_dividend(mf, static_cast<T>((modulus/2)-1));
}

template <typename M>
void test_mf_general_checks(const M& mf, typename M::IntegerType a,
                           typename M::IntegerType b, typename M::IntegerType c)
{
    namespace hc = ::hurchalla;
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

    T reference_two_a = tma::modadd(a, a, modulus);
    T reference_two_b = tma::modadd(b, b, modulus);
    EXPECT_TRUE(mf.convertOut(mf.two_times(x)) == reference_two_a);
    EXPECT_TRUE(mf.convertOut(mf.two_times(xc)) == reference_two_a);
    EXPECT_TRUE(mf.convertOut(mf.two_times(y)) == reference_two_b);
    EXPECT_TRUE(mf.convertOut(mf.two_times(yc)) == reference_two_b);
    EXPECT_TRUE(mf.getCanonicalValue(mf.two_times(x)) ==
                           mf.getCanonicalValue(mf.convertIn(reference_two_a)));
    EXPECT_TRUE(mf.getCanonicalValue(mf.two_times(xc)) ==
                           mf.getCanonicalValue(mf.convertIn(reference_two_a)));

    T diff1 = tma::modsub(b, a, modulus);
    test_subtract_variants(mf, y, x, diff1);
    T diff2 = tma::modsub(a, b, modulus);
    test_subtract_variants(mf, x, y, diff2);
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

    EXPECT_TRUE(mf.convertOut(mf.two_pow(0)) == 1);
    EXPECT_TRUE(mf.convertOut(mf.two_pow(1)) == 2);
    EXPECT_TRUE(mf.convertOut(mf.two_pow(3)) == tma::modpow<T>(2,3,modulus));
    EXPECT_TRUE(mf.convertOut(mf.two_pow(11)) == tma::modpow<T>(2,11,modulus));
    EXPECT_TRUE(mf.convertOut(mf.two_pow(127)) == tma::modpow<T>(2,127,modulus));

    EXPECT_TRUE(mf.convertOut(mf.pow(y,0)) == 1);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,1)) == b);
    EXPECT_TRUE(mf.convertOut(mf.pow(y,4)) == tma::modpow<T>(b,4,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,13)) == tma::modpow<T>(b,13,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,17)) == tma::modpow<T>(b,17,modulus));
    EXPECT_TRUE(mf.convertOut(mf.pow(y,127)) == tma::modpow<T>(b,127,modulus));

    // Do just a single test of pow()'s array form template function -
    // it's tested more thoroughly in test_montgomery_pow.cpp.
#ifdef __GNUC__
#  pragma GCC diagnostic push
// old versions of gcc and clang give unnecessary warnings about single braced
// initialization lists with std::array (newer versions fixed this).
#  pragma GCC diagnostic ignored "-Wmissing-braces"
#endif
    std::array<V,3> mv_bases = { x, y, z };
    std::array<T,3> t_bases = { a, b, c };
#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif
    // In case we're indirectly (or directly) using ConcreteMontgomeryForm,
    // ensure we use pow with an array size we've declared as possible.
    constexpr std::size_t possible_sizes[] = { TESTABLE_ARRAY_POW_SIZES() };
    constexpr std::size_t bases_size = std::tuple_size<decltype(mv_bases)>::value;
    // static_assert comparing to only possible_sizes[0] is a hack - to assert
    // properly, we'd write a constexpr function that searches possible_sizes.
    static_assert(possible_sizes[0] == bases_size, "");
    T exponent = 19;
    auto mv_res = mf.pow(mv_bases, exponent);
    for (std::size_t i = 0; i < bases_size; ++i) {
        auto correct_val = tma::modpow<T>(t_bases[i], exponent, modulus);
        EXPECT_TRUE(mf.convertOut(mv_res[i]) == correct_val);
    }
}




// The following is definitely not the typical way to create an instance of
// MontgomeryForm!  MontgomeryFactory exists solely to give our tests the
// option to use a run-time-polymorphic version of MontgomeryForm, to speed up
// compile times.  But the normal way to create an instance of MontgomeryForm
// is simply to call its constructor like you would any other class.
// An example of the normal way:
//    uint64_t modulus = 777777;
//    MontgomeryForm<uint64_t> mf(modulus);
//
template <class M, class ConcreteMF>
struct MontgomeryFactory {
    static M construct(typename M::IntegerType modulus)
    {
        auto* pc = new ConcreteMF(modulus);
        auto up = std::unique_ptr<typename ConcreteMF::Parent>(pc);
        return M(std::move(up));
    }
};
template <class M>
struct MontgomeryFactory<M, void> {
    static M construct(typename M::IntegerType modulus)
    {
        return M(modulus);
    }
};


template <typename M, class ConcreteMF = void>
void test_MontgomeryForm()
{
    namespace hc = ::hurchalla;
    using T = typename M::IntegerType;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;
    using FV = typename M::FusingValue;

    // see comments above MontgomeryFactory definition
    using MFactory = typename ::MontgomeryFactory<M, ConcreteMF>;

    // Try a basic test case first that is valid for all possible Monty types
    {
        T modulus = 13;
        // using MFactory is very unusual, but can help unit testing compile
        // times.  Normally we would have instead written
        //   M mf(modulus);
        M mf = MFactory::construct(modulus);

        V x = mf.convertIn(6);
        V y = mf.convertIn(11);

        C xc = mf.getCanonicalValue(x);
        C yc = mf.getCanonicalValue(y);

        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.two_times(x)) == 12);
        EXPECT_TRUE(mf.convertOut(mf.two_times(xc)) == 12);
        EXPECT_TRUE(mf.convertOut(mf.two_times(y)) == 9);
        EXPECT_TRUE(mf.convertOut(mf.two_times(yc)) == 9);
        test_subtract_variants(mf, y, x, 5);
        test_subtract_variants(mf, x, y, 8);
        T us = mf.convertOut(mf.unorderedSubtract(x,y));
        EXPECT_TRUE(us == 8 || us == 5);
        us = mf.convertOut(mf.unorderedSubtract(y,x));
        EXPECT_TRUE(us == 8 || us == 5);
        EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,y)) ==
                                         mf.getCanonicalValue(mf.convertIn(4)));
        EXPECT_TRUE(mf.getCanonicalValue(mf.add(x,yc)) ==
                                         mf.getCanonicalValue(mf.convertIn(4)));
        EXPECT_TRUE(mf.getCanonicalValue(mf.two_times(y)) ==
                                         mf.getCanonicalValue(mf.convertIn(9)));
        EXPECT_TRUE(mf.getCanonicalValue(mf.two_times(yc)) ==
                                         mf.getCanonicalValue(mf.convertIn(9)));
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
        M mf = MFactory::construct(modulus);
        V x = mf.convertIn(1);
        V y = mf.convertIn(2);

        C xc = mf.getCanonicalValue(x);
        C yc = mf.getCanonicalValue(y);

        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == 0);
        EXPECT_TRUE(mf.convertOut(mf.two_times(x)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.two_times(xc)) == 2);
        EXPECT_TRUE(mf.convertOut(mf.two_times(y)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.two_times(yc)) == 1);
        test_subtract_variants(mf, y, x, 1);
        test_subtract_variants(mf, x, y, 2);
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

    T max_modulus;
    {
        // normally we would just set max_modulus = M::max_modulus(),
        // but in the unusual case where M is an AbstractMontgomeryWrapper,
        // max_modulus() isn't static, and can't be static.
        M mf = MFactory::construct(static_cast<T>(3));  //we arbitrarily use 3 since it fits in any T
        max_modulus = mf.max_modulus();
    }

    // try the largest possible modulus
    {
        T modulus = max_modulus;
        M mf = MFactory::construct(modulus);
        EXPECT_TRUE(modulus > 4);

        V x = mf.convertIn(static_cast<T>(modulus - 1));
        V y = mf.convertIn(2);

        C xc = mf.getCanonicalValue(x);
        C yc = mf.getCanonicalValue(y);
    
        EXPECT_TRUE(mf.convertOut(mf.add(x,y)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.add(y,x)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.add(x,yc)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.add(y,xc)) == 1);
        EXPECT_TRUE(mf.convertOut(mf.two_times(x)) ==static_cast<T>(modulus-2));
        EXPECT_TRUE(mf.convertOut(mf.two_times(xc))==static_cast<T>(modulus-2));
        EXPECT_TRUE(mf.convertOut(mf.two_times(y)) == 4);
        EXPECT_TRUE(mf.convertOut(mf.two_times(yc)) == 4);
        EXPECT_TRUE(mf.two_times(xc) ==
                 mf.getCanonicalValue(mf.convertIn(static_cast<T>(modulus-2))));
        test_subtract_variants(mf, y, x, 3);
        test_subtract_variants(mf, x, y, modulus - 3);
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
        M mf = MFactory::construct(11);
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
        EXPECT_TRUE(max_modulus >= 5);
        M mf = MFactory::construct(static_cast<T>(max_modulus-2));
        T c = max_modulus-3;
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
        EXPECT_TRUE(max_modulus >= 5);
        M mf = MFactory::construct(static_cast<T>((max_modulus/4)*2 + 1));
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
        M mf = MFactory::construct(static_cast<T>(35));
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(28), GcdFunctor()) == 7);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(29), GcdFunctor()) == 1);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(70), GcdFunctor()) == 35);
    }
    if (117 <= max_modulus)
    {
        M mf = MFactory::construct(static_cast<T>(117));
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(78), GcdFunctor()) == 39);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(26), GcdFunctor()) == 13);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(27), GcdFunctor()) == 9);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(28), GcdFunctor()) == 1);
    }
    {
        M mf = MFactory::construct(static_cast<T>(3));  // smallest possible modulus
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(2), GcdFunctor()) == 1);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(0), GcdFunctor()) == 3);
    }
    {
        T modulus = max_modulus;
        ASSERT_TRUE(modulus > 9);
        while (modulus % 3 != 0 || modulus % 2 == 0)
            --modulus;
        M mf = MFactory::construct(modulus);
        EXPECT_TRUE(mf.gcd_with_modulus(mf.convertIn(12), GcdFunctor()) == 3);
    }

    // test remainder() and inverse() and divideBySmallPowerOf2()
    {
        T max = max_modulus;
        T mid = static_cast<T>(max/2);
        mid = (mid % 2 == 0) ? static_cast<T>(mid + 1) : mid;
        auto mf_3 = MFactory::construct(3);
        test_remainder(mf_3);    // smallest possible modulus
        test_inverse(mf_3);
        test_divideBySmallPowerOf2(mf_3);

        auto mf_max = MFactory::construct(max);
        test_remainder(mf_max);  // largest possible modulus
        test_inverse(mf_max);
        test_divideBySmallPowerOf2(mf_max);

        if (121 <= max) {
            auto mf_121 = MFactory::construct(121);
            test_remainder(mf_121);
            test_inverse(mf_121);
            test_divideBySmallPowerOf2(mf_121);
        }

        auto mf_mid = MFactory::construct(mid);
        test_remainder(mf_mid);
        test_inverse(mf_mid);
        test_divideBySmallPowerOf2(mf_mid);
    }
}

 
template <template<class,class> class MF, template<class> class MontyType>
void test_custom_monty()
{
    test_MontgomeryForm<MF<std::uint64_t, MontyType<std::uint64_t>>>();

#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
    test_MontgomeryForm<MF<std::uint8_t, MontyType<std::uint8_t>>>();
    test_MontgomeryForm<MF<std::uint16_t, MontyType<std::uint16_t>>>();
    test_MontgomeryForm<MF<std::uint32_t, MontyType<std::uint32_t>>>();
# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<MF<__uint128_t, MontyType<__uint128_t>>>();
# endif
#endif
}


} // end unnamed namespace

#endif

