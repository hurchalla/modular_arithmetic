// Copyright (c) 2020-2024 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

// Strictly for testing purposes, we'll ensure clockwork asserts are
// enabled by defining HURCHALLA_CLOCKWORK_ENABLE_ASSERTS.
#ifndef HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#  define HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#endif


#include "test_MontgomeryForm.h"

#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>


namespace {


namespace hc = ::hurchalla;

TEST(MontgomeryArithmetic, MontgomeryFormExamples) {

    // ---- Demonstrate simple modular addition ----
    {
        int64_t modulus = 15;
        int64_t x = 12;
        int64_t y = 4;

        // montgomery arithmetic always requires an odd modulus
        HPBC_CLOCKWORK_ASSERT(modulus % 2 == 1);
        // first construct a MontgomeryForm object to do Montgomery arithmetic
        // with a particular modulus.
        hc::MontgomeryForm<int64_t> mf(modulus);
        // convert x and y to their Montgomery representations.
        auto xm = mf.convertIn(x);
        auto ym = mf.convertIn(y);
        // perform modular addition on xm and ym
        auto summ = mf.add(xm, ym);
        // convert the Montgomery representation summ to normal integer domain
        int64_t result = mf.convertOut(summ);
        // Usually we keep values in the Montgomery domain for as long as we
        // can, and call convertOut only when we have finished all the modular
        // arithmetic that we wish to perform.  Note that this demonstration
        // works fine but is quite inefficient since it is doing just one
        // extremely simple operation in Montgomery domain - there is overhead
        // to construct the MontgomeryForm object and to call convertIn and
        // convertOut.  That is why you want to do many operations (or a CPU
        // intensive operation like modular exponentiation) in Montgomery
        // domain before calling convertOut; otherwise standard (non-Montgomery)
        // modular arithmetic will likely be more efficient.

        // Let's test that we get the same result using standard (non-
        // Montgomery) modular arithmetic.
        int64_t result2 = hc::modular_addition_prereduced_inputs(x, y, modulus);
        EXPECT_TRUE(result == result2);
    }

    // ---- Demonstrate modular exponentiation ----
    {
        int64_t modulus = 333333333;
        int64_t base = 42;
        int64_t exponent = 123456789;

        HPBC_CLOCKWORK_ASSERT(modulus % 2 == 1);  // montgomery arithmetic requires odd modulus
        // first construct a MontgomeryForm object to do Montgomery arithmetic
        // with a particular modulus.
        hc::MontgomeryForm<int64_t> mf(modulus);
        // convert base to its Montgomery representation.
        auto base_montval = mf.convertIn(base);
        // get the pow result in Montgomery representation.
        auto result_montval = mf.pow(base_montval, exponent);
        // convert the Montgomery representation result to normal integer domain
        int64_t result = mf.convertOut(result_montval);
        // Usually we keep values in the Montgomery domain for as long as we
        // can, and call convertOut only when we have finished all the modular
        // arithmetic that we wish to perform.

        // Let's test that we get the same result using standard (non-
        // Montgomery) modular arithmetic.  Note that Montgomery arithmetic is
        // usually faster whenever we have CPU intensive modular arithmetic,
        // which is why we use it.
        uint64_t result2 = hc::modular_pow(static_cast<uint64_t>(base),
                                           static_cast<uint64_t>(exponent),
                                           static_cast<uint64_t>(modulus));
        EXPECT_TRUE(static_cast<uint64_t>(result) == result2);
    }
}



// extensive tests of functionality with all possible Monty types ---



#if 1

// For unit testing, we want fast compile times, so it helps to use the version
// of MontgomeryForm that generally doesn't do force inlining.
#  if 1
constexpr bool forceInlineAllFunctions = false;
#  else
// note: even the default template arg for MontgomeryForm wouldn't have force
// inlined everything for uint128_t or int128_t (we would expect the functions
// to have too many instructions for it to be a good idea).  So for some T we
// get more inlining than the default, when this #else is enabled.
constexpr bool forceInlineAllFunctions = true;
#  endif
template <class T, class Monty> using MF = hc::MontgomeryForm<T, forceInlineAllFunctions, Monty>;
template <class T> using DefaultMF = hc::MontgomeryForm<T, forceInlineAllFunctions>;

#else
// this uses MontgomeryForm with the default amount of inlining.
template <class T, class Monty> using MF =
    hc::MontgomeryForm<T, (hurchalla::ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), Monty>;
template <class T> using DefaultMF =
    hc::MontgomeryForm<T>;
#endif



TEST(MontgomeryArithmetic, MontyQuarterRange) {
    test_custom_monty<MF, hc::detail::MontyQuarterRange>();
}

TEST(MontgomeryArithmetic, MontyHalfRange) {
    test_custom_monty<MF, hc::detail::MontyHalfRange>();
}


TEST(MontgomeryArithmetic, MontyFullRange) {
    namespace hcd = ::hurchalla::detail;
    test_MontgomeryForm<MF<uint8_t, hcd::MontyFullRange<std::uint8_t>>>();
    test_MontgomeryForm<MF<uint16_t, hcd::MontyFullRange<std::uint16_t>>>();
    test_MontgomeryForm<MF<uint32_t, hcd::MontyFullRange<std::uint32_t>>>();
    test_MontgomeryForm<MF<uint64_t, hcd::MontyFullRange<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<MF<__uint128_t, hcd::MontyFullRange<__uint128_t>>>();
#endif
}


TEST(MontgomeryArithmetic, MontyDefault) {
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

#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
// It would be absolutely normal and expected to use an unsigned integer type
// template argument for MontgomeryForm, but we can skip testing them here to
// save compilation time, because the resulting Montygomery classes using
// unsigned int types resolve to exactly the same types as will be tested in
// the MontyFullRange and MontyHalfRange and MontyQuarterRange TESTs above.
    test_MontgomeryForm<DefaultMF<std::uint8_t>>();
    test_MontgomeryForm<DefaultMF<std::uint16_t>>();
    test_MontgomeryForm<DefaultMF<std::uint32_t>>();
    test_MontgomeryForm<DefaultMF<std::uint64_t>>();
# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<DefaultMF<__uint128_t>>();
# endif
#endif

    test_MontgomeryForm<DefaultMF<std::int32_t>>();
#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
// To save compilation time we can also skip most signed integer type tests for
// plain MontgomeryForm.  These should differ from the unsigned versions (which
// in turn map to types we test above) only in the casts they perform for
// convertIn(), convertOut(), max_modulus(), getModulus() and gcd_with_modulus()
    test_MontgomeryForm<DefaultMF<std::int8_t>>();
    test_MontgomeryForm<DefaultMF<std::int16_t>>();
    test_MontgomeryForm<DefaultMF<std::int64_t>>();
// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<DefaultMF<__int128_t>>();
# endif
#endif
}


} // end anonymous namespace
