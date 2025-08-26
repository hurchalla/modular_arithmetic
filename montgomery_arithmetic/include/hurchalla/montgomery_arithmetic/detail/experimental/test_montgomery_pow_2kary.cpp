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

#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
// You can define the next macro to test all realistic function template arg
// combos, rather than only the template args likely to be used.
// This macro can result in very long compile and execution times...
//#define HURCHALLA_TEST_POW_2KARY_ULTRA_HEAVYWEIGHT 1
#endif


#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/detail/impl_montgomery_pow_2kary.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>
#include <array>
#include <utility>
#include <vector>
#include <iostream>
#include <random>

#include<string>


namespace {


namespace hc = ::hurchalla;
using pow_2kary = ::hurchalla::detail::impl_montgomery_pow_2kary;


/*
template <typename U>
std::string uint_to_string(U number)
{
   namespace hc = hurchalla;
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   if (number == 0)
      return std::string("0");
   std::string str;
   while (number > 0) {
      char digit = static_cast<char>((number % 10) + '0');
      str.push_back(digit);
      number = number / 10;
   }
   return std::string(str.rbegin(), str.rend());
}

    std::cout << "modulus = " << uint_to_string(modulus) << ", exponent = " << uint_to_string(exponent)
              << ", answer = " << uint_to_string(answer) << ", result = " << uint_to_string(result) << "\n";
*/



template <std::size_t NUM_BASES, class M, typename U>
void test_pow_2kary_array(M& mf, typename M::IntegerType base, U exponent)
{
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");

    using T = typename M::IntegerType;
    using V = typename M::MontgomeryValue;

    T modulus = mf.getModulus();

    std::array<V, NUM_BASES> mv_bases;
    std::array<T, NUM_BASES> answers;
    for (std::size_t i=0; i<mv_bases.size(); ++i) {
        T tmpbase = static_cast<T>((base + i) % modulus);
        mv_bases[i] = mf.convertIn(tmpbase);
        answers[i] = hc::modular_pow(tmpbase, exponent, modulus);
    }

    std::array<V, NUM_BASES> results;
#ifndef HURCHALLA_TEST_POW_2KARY_ULTRA_HEAVYWEIGHT
    results = pow_2kary::call<M,U,NUM_BASES,true,4>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,false,5>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
#else
    results = pow_2kary::call<M,U,NUM_BASES,true,1>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,true,2>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,true,3>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,true,4>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,true,5>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,true,6>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,true,7>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }

    results = pow_2kary::call<M,U,NUM_BASES,false,1>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,false,2>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,false,3>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,false,4>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,false,5>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,false,6>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
    results = pow_2kary::call<M,U,NUM_BASES,false,7>(mf, mv_bases, exponent);
    for (std::size_t i=0; i<results.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(results[i]) == answers[i]);
    }
#endif
}



// this utility function vector_to_array() is adapted from
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
std::array<T, sizeof...(N)> vector_to_array_impl(std::vector<T>& vec, seq<N...>)
{
   HPBC_CLOCKWORK_PRECONDITION2(vec.size() >= sizeof...(N));
   return { vec[N]... };
}
template<std::size_t SIZE, class T>
std::array<T, SIZE> vector_to_array(std::vector<T>& vec)
{
   HPBC_CLOCKWORK_PRECONDITION2(vec.size() >= SIZE);
   return vector_to_array_impl(vec, typename genseq<SIZE>::type{} );
}



template <class MF, std::size_t ARRAY_SIZE, typename U>
void test_pow_2kary_allarrays(typename MF::IntegerType starting_modulus,
                    typename MF::IntegerType starting_base, U starting_exponent)
{
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
    using T = typename MF::IntegerType;
    using V = typename MF::MontgomeryValue;

    using ULL = unsigned long long;
    static std::mt19937_64 gen;
    std::uniform_int_distribution<ULL> distribU;

    std::vector<MF> mf_vec;
    std::array<T, ARRAY_SIZE> bases;
    std::array<V, ARRAY_SIZE> mv_bases;
    std::array<U, ARRAY_SIZE> exponents;
    T m = starting_modulus;
    T base = starting_base;
    U expo = starting_exponent;
    for (std::size_t i=0; i<ARRAY_SIZE; ++i) {
        mf_vec.push_back(MF(m));
        bases[i] = base % m;
        mv_bases[i] = mf_vec[i].convertIn(bases[i]);
        exponents[i] = expo;

        if (m < 3)
            m = 3;
        else if (m >= MF::max_modulus() - 1)
            m = MF::max_modulus();
        else
            m += 2;
        ULL randnum = distribU(gen);
        if (i % 2 == 0)
            expo = static_cast<U>(expo + randnum);
        else
            expo = static_cast<U>(expo - randnum);
        randnum = distribU(gen);
        if (i % 2 == 0)
            base = static_cast<T>(base + randnum);
        else
            base = static_cast<T>(base - randnum);
    }
    std::array<MF, ARRAY_SIZE> mfs = vector_to_array<ARRAY_SIZE>(mf_vec);

    std::array<T, ARRAY_SIZE> answers;
    for (std::size_t i=0; i<ARRAY_SIZE; ++i) {
        answers[i] = hc::modular_pow<T>(
                                   bases[i], exponents[i], mfs[i].getModulus());
    }

    std::array<V, ARRAY_SIZE> mv_results;
#ifndef HURCHALLA_TEST_POW_2KARY_ULTRA_HEAVYWEIGHT
    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,4>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);
#else
    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,1>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);

    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,2>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);

    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,3>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);

    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,4>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);

    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,5>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);

    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,6>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);

    mv_results = pow_2kary::call<MF,U,ARRAY_SIZE,7>(mfs, mv_bases, exponents);
    for (std::size_t i=0; i<ARRAY_SIZE; ++i)
        EXPECT_TRUE(mfs[i].convertOut(mv_results[i]) == answers[i]);
#endif
}




template <typename M, typename U>
void test_pow_2kary(typename M::IntegerType modulus,
                    typename M::IntegerType base, U exponent)
{
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
    using T = typename M::IntegerType;
    using V = typename M::MontgomeryValue;

    M mf(modulus);

    // first try the non-array overload of pow
    T answer = hc::modular_pow<T>(base, exponent, modulus);
    V mv_base = mf.convertIn(base);
    V mv_result;
#ifndef HURCHALLA_TEST_POW_2KARY_ULTRA_HEAVYWEIGHT
    mv_result = pow_2kary::call<M, U, false, 4>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, true, 5>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
#else
    mv_result = pow_2kary::call<M, U, false, 1>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, false, 2>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, false, 3>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, false, 4>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, false, 5>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, false, 6>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, false, 7>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);

    mv_result = pow_2kary::call<M, U, true, 1>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, true, 2>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, true, 3>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, true, 4>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, true, 5>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, true, 6>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
    mv_result = pow_2kary::call<M, U, true, 7>(mf, mv_base, exponent);
    EXPECT_TRUE(mf.convertOut(mv_result) == answer);
#endif

    // test the (partial) array version of pow_2kary using different array sizes
    test_pow_2kary_array<2>(mf, base, exponent);
#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
    test_pow_2kary_array<1>(mf, base, exponent);
    test_pow_2kary_array<3>(mf, base, exponent);
    test_pow_2kary_array<4>(mf, base, exponent);
    test_pow_2kary_array<9>(mf, base, exponent);
    //test_pow_2kary_array<14>(mf, base, exponent);
    //test_pow_2kary_array<29>(mf, base, exponent);
    //test_pow_2kary_array<61>(mf, base, exponent);
    //test_pow_2kary_array<120>(mf, base, exponent);
#endif

    // test the all-array param version of pow_2kary using different array sizes
    test_pow_2kary_allarrays<M,2>(modulus, base, exponent);
#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
    test_pow_2kary_allarrays<M,1>(modulus, base, exponent);
    test_pow_2kary_allarrays<M,3>(modulus, base, exponent);
    test_pow_2kary_allarrays<M,5>(modulus, base, exponent);
    test_pow_2kary_allarrays<M,8>(modulus, base, exponent);
    //test_pow_2kary_allarrays<M,13>(modulus, base, exponent);
    //test_pow_2kary_allarrays<M,23>(modulus, base, exponent);
    //test_pow_2kary_allarrays<M,51>(modulus, base, exponent);
    //test_pow_2kary_allarrays<M,112>(modulus, base, exponent);
#endif
}



template <typename M, typename U>
void run_pow_tests()
{
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
    using T = typename M::IntegerType;

    // Try a basic test case first that is valid for all possible Monty types
    {
        T modulus = 15;
        T base = 8;
        U exponent = 17;
        test_pow_2kary<M>(modulus, base, exponent);
    }
    // Try a test with the smallest possible modulus
    {
        T modulus = 3;
        T base = 2;
        U exponent = 9;
        test_pow_2kary<M>(modulus, base, exponent);
    }
    // Try the largest possible modulus
    {
        T modulus = M::max_modulus();
        T base = static_cast<T>(modulus - 1);
        U exponent = 183;
        test_pow_2kary<M>(modulus, base, exponent);
    }


    // Try a bunch of general tests...

    if (M::max_modulus() >= 119) {
        T modulus = 119;
        T base = 5;
        U exponent = 6;
        test_pow_2kary<M>(modulus, base, exponent);
        base = 10;
        exponent = 0;
        test_pow_2kary<M>(modulus, base, exponent);
        base = 0;
        exponent = 0;
        test_pow_2kary<M>(modulus, base, exponent);
        base = modulus - 1;
        exponent = 1;
        test_pow_2kary<M>(modulus, base, exponent);
#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4309)
#endif
        base = 0;
        exponent = static_cast<U>(1326);
        test_pow_2kary<M>(modulus, base, exponent);
        base = 1;
        exponent = static_cast<U>(551);
        test_pow_2kary<M>(modulus, base, exponent);
        base = 65;
        exponent = 1;
        test_pow_2kary<M>(modulus, base, exponent);
        base = 73;
        exponent = static_cast<U>(933);
        test_pow_2kary<M>(modulus, base, exponent);
#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
    }
    {
        T max = M::max_modulus();
        T modulus = static_cast<T>(max-2);
        T base = static_cast<T>(max-3);
        U exponent = 22;
        test_pow_2kary<M>(modulus, base, exponent);
        base = static_cast<T>(max/2);
        exponent = 49;
        test_pow_2kary<M>(modulus, base, exponent);
        base = static_cast<T>(max/2-1);
        exponent = 252;
        test_pow_2kary<M>(modulus, base, exponent);
        exponent = hc::ut_numeric_limits<U>::max();
        test_pow_2kary<M>(modulus, base, exponent);
        base = 1;
        exponent = 125;
        test_pow_2kary<M>(modulus, base, exponent);
    }
    {
        T modulus = (M::max_modulus()/4)*2 + 1;
        T base = 5;
        U exponent = 89;
        test_pow_2kary<M>(modulus, base, exponent);
        base = static_cast<T>(modulus-4);
        exponent = 3;
        test_pow_2kary<M>(modulus, base, exponent);
        base = static_cast<T>(modulus/2);
        exponent = 2;
        test_pow_2kary<M>(modulus, base, exponent);
        base = static_cast<T>(modulus/2-1);
        exponent = 4;
        test_pow_2kary<M>(modulus, base, exponent);
        base = 0;
        exponent = 123;
        test_pow_2kary<M>(modulus, base, exponent);
        base = 7;
        exponent = hc::ut_numeric_limits<U>::max();
        test_pow_2kary<M>(modulus, base, exponent);
        base = static_cast<T>(modulus - 5);
        --exponent;
        test_pow_2kary<M>(modulus, base, exponent);
        exponent = exponent/2;
        test_pow_2kary<M>(modulus, base, exponent);
    }
}



// For unit testing, we want fast compile times, so it helps to use the version
// of MontgomeryForm that generally doesn't do force inlining.
#if 1
constexpr bool forceInlineAllFunctions = false;
#else
// note: even the default template arg for MontgomeryForm wouldn't have force
// inlined everything for uint128_t or int128_t (we would expect the functions
// to have too many instructions for it to be a good idea).  So for some T we
// get more inlining than the default, when this #else is enabled.
constexpr bool forceInlineAllFunctions = true;
#endif

template <class T, class Monty> using MF = hc::MontgomeryForm<T, forceInlineAllFunctions, Monty>;
template <class T> using DefaultMF = hc::MontgomeryForm<T, forceInlineAllFunctions>;



TEST(MontgomeryArithmetic, montgomery_pow_2kary) {

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
    run_pow_tests<MF<std::uint64_t,
                    hc::detail::MontyWrappedStandardMath<std::uint64_t>>, U2>();

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
