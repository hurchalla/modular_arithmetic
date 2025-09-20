// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// Strictly for testing purposes,  we'll define HURCHALLA_ALLOW_INLINE_ASM_ALL,
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


#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/detail/impl_montgomery_two_pow.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>
#include <array>
#include <utility>
#include <vector>
#include <random>

namespace {


namespace hc = ::hurchalla;
using two_pow = ::hurchalla::detail::impl_montgomery_two_pow;


// this utility function vector_to_stdarray() is adapted from
// https://stackoverflow.com/a/18497366
template<std::size_t ... N>
struct seq
{
   using type = seq<N...>;
   template <std::size_t I> struct push_back : seq<N..., I> {};
};
template<std::size_t N>
struct genseq : genseq<N-1>::type::template push_back<N-1> {};
template<>
struct genseq<0> : seq<> {};

template<class T, std::size_t... N>
std::array<T, sizeof...(N)> vector_to_stdarray_impl(const std::vector<T>& vec, seq<N...>)
{
   HPBC_CLOCKWORK_PRECONDITION(vec.size() >= sizeof...(N));
   return { vec[N]... };
}
template<std::size_t SIZE, class T>
std::array<T, SIZE> vector_to_stdarray(const std::vector<T>& vec)
{
   HPBC_CLOCKWORK_PRECONDITION(vec.size() >= SIZE);
   return vector_to_stdarray_impl(vec, typename genseq<SIZE>::type{} );
}



template <class MF, std::size_t ARRAY_SIZE, typename U>
void test_two_pow_array(typename MF::IntegerType starting_modulus,
                        U starting_exponent)
{
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
    using T = typename MF::IntegerType;
    using V = typename MF::MontgomeryValue;

    std::vector<MF> mf_vec;
    std::array<U, ARRAY_SIZE> exponents;
    T m = starting_modulus;
    U expo = starting_exponent;
    for (std::size_t i=0; i<ARRAY_SIZE; ++i) {
        mf_vec.push_back(MF(m));
        exponents[i] = expo;
        if (m < 3)
            m = 3;
        else if (m >= MF::max_modulus() - 1)
            m = MF::max_modulus();
        else
            m += 2;
        ++expo;
    }
    std::array<MF, ARRAY_SIZE> mfs = vector_to_stdarray<ARRAY_SIZE>(mf_vec);

    std::array<V,ARRAY_SIZE> results;

    results = two_pow::call<0,0,MF,U,ARRAY_SIZE>(mfs, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i) {
        EXPECT_TRUE(mfs[i].convertOut(results[i]) ==
                      hc::modular_pow<T>(2, exponents[i], mfs[i].getModulus()));
    }
    results = two_pow::call<0,2,MF,U,ARRAY_SIZE>(mfs, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i) {
        EXPECT_TRUE(mfs[i].convertOut(results[i]) ==
                      hc::modular_pow<T>(2, exponents[i], mfs[i].getModulus()));
    }
}



template <typename M, typename U>
void test_two_pow(typename M::IntegerType modulus, U exponent)
{
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
    using T = typename M::IntegerType;

    // first test the non-array two_pow
    M mf(modulus);
    T answer = hc::modular_pow<T>(2, exponent, modulus);
    T result;
    result = mf.convertOut(two_pow::call<true,0,2,M,U>(mf,exponent));
    EXPECT_TRUE(result == answer);
    result = mf.convertOut(two_pow::call<true,0,3,M,U>(mf,exponent));
    EXPECT_TRUE(result == answer);
    result = mf.convertOut(two_pow::call<false,0,2,M,U>(mf,exponent));
    EXPECT_TRUE(result == answer);
    result = mf.convertOut(two_pow::call<false,0,3,M,U>(mf,exponent));
    EXPECT_TRUE(result == answer);

    // test the array version of two_pow with different array sizes
    test_two_pow_array<M,1>(modulus, exponent);
    test_two_pow_array<M,2>(modulus, exponent);
    test_two_pow_array<M,3>(modulus, exponent);
#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
    test_two_pow_array<M,4>(modulus, exponent);
    test_two_pow_array<M,5>(modulus, exponent);
    test_two_pow_array<M,6>(modulus, exponent);
    test_two_pow_array<M,7>(modulus, exponent);
    test_two_pow_array<M,8>(modulus, exponent);

    test_two_pow_array<M,9>(modulus, exponent);
    test_two_pow_array<M,10>(modulus, exponent);
    test_two_pow_array<M,11>(modulus, exponent);
    test_two_pow_array<M,12>(modulus, exponent);
    test_two_pow_array<M,13>(modulus, exponent);
    test_two_pow_array<M,19>(modulus, exponent);
    //test_two_pow_array<M,61>(modulus, exponent);
    //test_two_pow_array<M,120>(modulus, exponent);
#endif
}



template <typename U>
U generate_random_value(std::mt19937_64& gen,
                        std::uniform_int_distribution<uint64_t>& distrib64)
{
   static_assert(hurchalla::ut_numeric_limits<U>::is_integer, "");
   static_assert(!hurchalla::ut_numeric_limits<U>::is_signed, "");
   static_assert(hurchalla::ut_numeric_limits<U>::digits <= 128, "");
   if HURCHALLA_CPP17_CONSTEXPR (hurchalla::ut_numeric_limits<U>::digits > 64) {
      uint64_t u1 = distrib64(gen);
      uint64_t u2 = distrib64(gen);
      using P = typename hurchalla::safely_promote_unsigned<U>::type;
      U val = static_cast<U>((static_cast<P>(u2) << 64u) | u1);
      return val;
   } else {
      return static_cast<U>(distrib64(gen));
   }
}


template <typename M, typename U>
void run_pow_tests()
{
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
    using T = typename M::IntegerType;

    std::mt19937_64 gen(2);   // 2 is an arbitrary seed
    std::uniform_int_distribution<uint64_t> distrib64;

    // Try a basic test case first that is valid for all possible Monty types
    {
        T modulus = 13;
        U exponent = 11;
        test_two_pow<M>(modulus, exponent);
    }
    // Try a test with the smallest possible modulus
    {
        T modulus = 3;
        U exponent = 5;
        test_two_pow<M>(modulus, exponent);
    }
    // Try the largest possible modulus
    {
        T modulus = M::max_modulus();
        U exponent = 179;
        test_two_pow<M>(modulus, exponent);
    }


    // Try a bunch of general tests...

    if (M::max_modulus() >= 113) {
        T modulus = 113;
        U exponent = 6;
        test_two_pow<M>(modulus, exponent);
        exponent = 0;
        test_two_pow<M>(modulus, exponent);
        exponent = 1;
        test_two_pow<M>(modulus, exponent);
        exponent = 7;
        test_two_pow<M>(modulus, exponent);
        exponent = 8;
        test_two_pow<M>(modulus, exponent);
#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4309)
#endif
        exponent = static_cast<U>(1356);
        test_two_pow<M>(modulus, exponent);
        exponent = static_cast<U>(541);
        test_two_pow<M>(modulus, exponent);
        exponent = static_cast<U>(934);
        test_two_pow<M>(modulus, exponent);
#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
    }
    {
        T max = M::max_modulus();
        T modulus = static_cast<T>(max-2);
        U exponent = 24;
        test_two_pow<M>(modulus, exponent);
        exponent = 43;
        test_two_pow<M>(modulus, exponent);
        exponent = 253;
        test_two_pow<M>(modulus, exponent);
        exponent = 135;
        test_two_pow<M>(modulus, exponent);
        exponent = hc::ut_numeric_limits<U>::max();
        test_two_pow<M>(modulus, exponent);
        exponent -= 4;
        test_two_pow<M>(modulus, exponent);

        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
    }
    {
        T modulus = (M::max_modulus()/4)*2 + 1;
        U exponent = 89;
        test_two_pow<M>(modulus, exponent);
        exponent = 3;
        test_two_pow<M>(modulus, exponent);
        exponent = 2;
        test_two_pow<M>(modulus, exponent);
        exponent = 123;
        test_two_pow<M>(modulus, exponent);
        exponent = hc::ut_numeric_limits<U>::max();
        test_two_pow<M>(modulus, exponent);
        --exponent;
        test_two_pow<M>(modulus, exponent);

        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
    }
    {
        using TU = typename hc::extensible_make_unsigned<T>::type;
        T modulus = static_cast<T>(generate_random_value<TU>(gen, distrib64));
        while (modulus <= 10)
            modulus = static_cast<T>(generate_random_value<TU>(gen, distrib64));
        while (modulus > M::max_modulus())
            modulus /= 2;
        modulus = static_cast<T>(modulus + (modulus % 2) - 1);

        U exponent = hc::ut_numeric_limits<U>::max();
        test_two_pow<M>(modulus, exponent);
        exponent -= 123;
        test_two_pow<M>(modulus, exponent);

        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
    }
    {
        using TU = typename hc::extensible_make_unsigned<T>::type;
        T modulus = static_cast<T>(generate_random_value<TU>(gen, distrib64));
        while (modulus <= 10)
            modulus = static_cast<T>(generate_random_value<TU>(gen, distrib64));
        while (modulus > M::max_modulus())
            modulus /= 2;
        modulus = static_cast<T>(modulus + (modulus % 2) - 1);

        U exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
        exponent = generate_random_value<U>(gen, distrib64);
        test_two_pow<M>(modulus, exponent);
    }
}



// For unit testing, we want fast compile times, so it helps to use the version
// of MontgomeryForm that doesn't do force inlining.
#if 1
constexpr bool forceInlineAllFunctions = false;
#else
constexpr bool forceInlineAllFunctions = true;
#endif

template <class T, class Monty> using MF = hc::MontgomeryForm<T, forceInlineAllFunctions, Monty>;
template <class T> using DefaultMF = hc::MontgomeryForm<T, forceInlineAllFunctions>;



TEST(MontgomeryArithmetic, montgomery_two_pow) {

    run_pow_tests<DefaultMF<std::uint8_t>, std::uint16_t>();

    using U1 = std::uint16_t;
    run_pow_tests<MF<std::uint8_t,
                    hc::detail::MontyQuarterRange<std::uint8_t>>, U1>();
    run_pow_tests<MF<std::uint8_t,
                    hc::detail::MontyHalfRange<std::uint8_t>>, U1>();
    run_pow_tests<MF<std::uint8_t,
                    hc::detail::MontyFullRange<std::uint8_t>>, U1>();
    run_pow_tests<MF<std::uint8_t,
                    hc::detail::MontyWrappedStandardMath<std::uint8_t>>, U1>();

    run_pow_tests<MF<std::uint16_t,
                    hc::detail::MontyQuarterRange<std::uint16_t>>, U1>();
    run_pow_tests<MF<std::uint16_t,
                    hc::detail::MontyHalfRange<std::uint16_t>>, U1>();
    run_pow_tests<MF<std::uint16_t,
                    hc::detail::MontyFullRange<std::uint16_t>>, U1>();
    run_pow_tests<MF<std::uint16_t,
                    hc::detail::MontyWrappedStandardMath<std::uint16_t>>, U1>();

    using U2 = std::uint64_t;
    run_pow_tests<MF<std::uint32_t,
                    hc::detail::MontyQuarterRange<std::uint32_t>>, U2>();
    run_pow_tests<MF<std::uint32_t,
                    hc::detail::MontyHalfRange<std::uint32_t>>, U2>();
    run_pow_tests<MF<std::uint32_t,
                    hc::detail::MontyFullRange<std::uint32_t>>, U2>();
    run_pow_tests<MF<std::uint32_t,
                    hc::detail::MontyWrappedStandardMath<std::uint32_t>>, U2>();

    run_pow_tests<MF<std::uint64_t,
                    hc::detail::MontyQuarterRange<std::uint64_t>>, U2>();
    run_pow_tests<MF<std::uint64_t,
                    hc::detail::MontyHalfRange<std::uint64_t>>, U2>();
    run_pow_tests<MF<std::uint64_t,
                    hc::detail::MontyFullRange<std::uint64_t>>, U2>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    using U3 = __uint128_t;
#else
    using U3 = std::uint64_t;
#endif
    run_pow_tests<MF<std::uint64_t,
                    hc::detail::MontyWrappedStandardMath<std::uint64_t>>, U3>();

#if HURCHALLA_COMPILER_HAS_UINT128_T()
    using U4 = __uint128_t;
    run_pow_tests<MF<__uint128_t,
                    hc::detail::MontyQuarterRange<__uint128_t>>, U4>();
    run_pow_tests<MF<__uint128_t,
                    hc::detail::MontyHalfRange<__uint128_t>>, U4>();
    run_pow_tests<MF<__uint128_t,
                    hc::detail::MontyFullRange<__uint128_t>>, U4>();
    run_pow_tests<MF<__uint128_t,
                    hc::detail::MontyWrappedStandardMath<__uint128_t>>, U4>();
#endif
}


} // end unnamed namespace
