// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontgomeryFormExtensions.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontyFullRangeMasked.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/ConcreteMontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryWrapper.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>
#include <cassert>
#include <random>


namespace {


template <typename U>
U generate_random_value(std::mt19937_64& gen,
                        std::uniform_int_distribution<uint64_t>& distrib64)
{
   static_assert(hurchalla::ut_numeric_limits<U>::is_integer, "");
   static_assert(!hurchalla::ut_numeric_limits<U>::is_signed, "");
   static_assert(hurchalla::ut_numeric_limits<U>::digits <= 128, "");
   if (hurchalla::ut_numeric_limits<U>::digits > 64) {
      uint64_t u1 = distrib64(gen);
      uint64_t u2 = distrib64(gen);
      using P = typename hurchalla::safely_promote_unsigned<U>::type;
      U val = static_cast<U>((static_cast<P>(u2) << 64u) | u1);
      return val;
   } else {
      return static_cast<U>(distrib64(gen));
   }
}



template <class MF>
void MFE_tests(const MF& mf)
{
    namespace hc = ::hurchalla;
    using T = typename MF::IntegerType;
    using V = typename MF::MontgomeryValue;
    using C = typename MF::CanonicalValue;

    using MFE0 = hc::detail::MontgomeryFormExtensions<MF, hc::LowlatencyTag>;
    using MFE1 = hc::detail::MontgomeryFormExtensions<MF, hc::LowuopsTag>;

    using RU = typename MFE1::RU;
    static_assert(hc::ut_numeric_limits<RU>::is_integer, "");
    static_assert(!hc::ut_numeric_limits<RU>::is_signed, "");

    constexpr RU maxRU = hc::ut_numeric_limits<RU>::max();
    constexpr int digitsR = hc::ut_numeric_limits<RU>::digits;

    C mont_one = mf.getUnityValue();
    C mont_two = mf.add(mont_one, mont_one);
    V mont_sqrtR = mf.pow(mf.convertIn(2), static_cast<T>(digitsR/2));
    V mont_R = mf.square(mont_sqrtR);


    // test getMagicValue
    // (magic should be R cubed mod N  (in normal integer form))
    RU magic0 = MFE0::getMagicValue(mf);
    RU magic1 = MFE1::getMagicValue(mf);

    // test convertInExtended
    {
        C c0, c1, cAnswer;
        RU a = 3;
        c0 = mf.getCanonicalValue(MFE0::convertInExtended(mf, a));
        c1 = mf.getCanonicalValue(MFE1::convertInExtended(mf, a));
        cAnswer = mf.getCanonicalValue(mf.convertIn(static_cast<T>(a)));
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);

        a = maxRU;
        c0 = mf.getCanonicalValue(MFE0::convertInExtended(mf, a));
        c1 = mf.getCanonicalValue(MFE1::convertInExtended(mf, a));
        cAnswer = mf.getCanonicalValue(mf.subtract(mont_R, mont_one));
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);
    }
    // test convertInExtended_aTimesR
    {
        C c0, c1, cAnswer;
        RU a = 7;
        c0 = mf.getCanonicalValue(MFE0::convertInExtended_aTimesR(mf,a,magic0));
        c1 = mf.getCanonicalValue(MFE1::convertInExtended_aTimesR(mf,a,magic1));
        V mont_a = mf.convertIn(static_cast<T>(a));
        cAnswer = mf.getCanonicalValue(mf.multiply(mont_a, mont_R));
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);

        a = maxRU - 1;
        c0 = mf.getCanonicalValue(MFE0::convertInExtended_aTimesR(mf,a,magic0));
        c1 = mf.getCanonicalValue(MFE1::convertInExtended_aTimesR(mf,a,magic1));
        mont_a = mf.subtract(mont_R, mont_two);
        cAnswer = mf.getCanonicalValue(mf.multiply(mont_a, mont_R));
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);
    }
    // test twoPowLimited
    {
        C c0, c1, cAnswer;
        size_t exponent;
        exponent = 0;
        c0 = mf.getCanonicalValue(MFE0::twoPowLimited(mf, exponent));
        c1 = mf.getCanonicalValue(MFE1::twoPowLimited(mf, exponent));
        V tmp = mf.pow(mf.convertIn(2), static_cast<T>(exponent));
        cAnswer = mf.getCanonicalValue(tmp);
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);

        exponent = 7;
        c0 = mf.getCanonicalValue(MFE0::twoPowLimited(mf, exponent));
        c1 = mf.getCanonicalValue(MFE1::twoPowLimited(mf, exponent));
        tmp = mf.pow(mf.convertIn(2), static_cast<T>(exponent));
        cAnswer = mf.getCanonicalValue(tmp);
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);

        exponent = static_cast<size_t>(digitsR - 5);
        c0 = mf.getCanonicalValue(MFE0::twoPowLimited(mf, exponent));
        c1 = mf.getCanonicalValue(MFE1::twoPowLimited(mf, exponent));
        tmp = mf.pow(mf.convertIn(2), static_cast<T>(exponent));
        cAnswer = mf.getCanonicalValue(tmp);
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);

        exponent = static_cast<size_t>(digitsR - 1);
        c0 = mf.getCanonicalValue(MFE0::twoPowLimited(mf, exponent));
        c1 = mf.getCanonicalValue(MFE1::twoPowLimited(mf, exponent));
        tmp = mf.pow(mf.convertIn(2), static_cast<T>(exponent));
        cAnswer = mf.getCanonicalValue(tmp);
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);
    }
    // test RTimesTwoPowLimited
    {
        C c0, c1, cAnswer;
        size_t expnt;   // exponent
        expnt = 4;
        c0 = mf.getCanonicalValue(MFE0::RTimesTwoPowLimited(mf, expnt, magic0));
        c1 = mf.getCanonicalValue(MFE1::RTimesTwoPowLimited(mf, expnt, magic1));
        V tmp = mf.pow(mf.convertIn(2), static_cast<T>(expnt));
        cAnswer = mf.getCanonicalValue(mf.multiply(tmp, mont_R));
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);

        expnt = static_cast<size_t>(digitsR - 2);
        c0 = mf.getCanonicalValue(MFE0::RTimesTwoPowLimited(mf, expnt, magic0));
        c1 = mf.getCanonicalValue(MFE1::RTimesTwoPowLimited(mf, expnt, magic1));
        tmp = mf.pow(mf.convertIn(2), static_cast<T>(expnt));
        cAnswer = mf.getCanonicalValue(mf.multiply(tmp, mont_R));
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);

        expnt = static_cast<size_t>(digitsR - 1);
        c0 = mf.getCanonicalValue(MFE0::RTimesTwoPowLimited(mf, expnt, magic0));
        c1 = mf.getCanonicalValue(MFE1::RTimesTwoPowLimited(mf, expnt, magic1));
        tmp = mf.pow(mf.convertIn(2), static_cast<T>(expnt));
        cAnswer = mf.getCanonicalValue(mf.multiply(tmp, mont_R));
        EXPECT_TRUE(c0 == cAnswer);
        EXPECT_TRUE(c1 == cAnswer);
    }
}





template <typename M>
void test_MFE()
{
    namespace hc = ::hurchalla;
    using T = typename M::IntegerType;

    std::mt19937_64 gen(2);   // 2 is an arbitrary seed
    std::uniform_int_distribution<uint64_t> distrib64;

    using MFE = hc::detail::MontgomeryFormExtensions<M, hc::LowlatencyTag>;
    using RU = typename MFE::RU;

    // Try a basic test case first that is valid for all possible Monty types
    {
        T modulus = 11;
        M mf(modulus);
        MFE_tests(mf);
    }

    T max_modulus = M::max_modulus();

    {
        T modulus = max_modulus;
        M mf(modulus);
        MFE_tests(mf);
    }
    {
        T modulus = max_modulus - 18;
        M mf(modulus);
        MFE_tests(mf);
    }
    {
        T modulus = static_cast<T>(generate_random_value<RU>(gen, distrib64));
        while (modulus < 10)
            modulus = static_cast<T>(generate_random_value<RU>(gen, distrib64));
        while (modulus > max_modulus)
            modulus /= 2;
        modulus = static_cast<T>(modulus + (modulus % 2) - 1);
        M mf(modulus);
        MFE_tests(mf);
    }
    {
        T modulus = static_cast<T>(generate_random_value<RU>(gen, distrib64));
        while (modulus < 10)
            modulus = static_cast<T>(generate_random_value<RU>(gen, distrib64));
        while (modulus > max_modulus)
            modulus /= 2;
        modulus = static_cast<T>(modulus + (modulus % 2) - 1);
        M mf(modulus);
        MFE_tests(mf);
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

template <class T, class Monty> using MF =
    hurchalla::MontgomeryForm<T, forceInlineAllFunctions, Monty>;



TEST(MontgomeryFormExtensions, MontyQuarterRange) {
    namespace hcd = ::hurchalla::detail;
    test_MFE<MF<uint8_t, hcd::MontyQuarterRange<std::uint8_t>>>();
    test_MFE<MF<uint16_t, hcd::MontyQuarterRange<std::uint16_t>>>();
    test_MFE<MF<uint32_t, hcd::MontyQuarterRange<std::uint32_t>>>();
    test_MFE<MF<uint64_t, hcd::MontyQuarterRange<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MFE<MF<__uint128_t, hcd::MontyQuarterRange<__uint128_t>>>();
#endif
}

TEST(MontgomeryFormExtensions, MontyHalfRange) {
    namespace hcd = ::hurchalla::detail;
    test_MFE<MF<uint8_t, hcd::MontyHalfRange<std::uint8_t>>>();
    test_MFE<MF<uint16_t, hcd::MontyHalfRange<std::uint16_t>>>();
    test_MFE<MF<uint32_t, hcd::MontyHalfRange<std::uint32_t>>>();
    test_MFE<MF<uint64_t, hcd::MontyHalfRange<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MFE<MF<__uint128_t, hcd::MontyHalfRange<__uint128_t>>>();
#endif
}

TEST(MontgomeryFormExtensions, MontyFullRange) {
    namespace hcd = ::hurchalla::detail;
    test_MFE<MF<uint8_t, hcd::MontyFullRange<std::uint8_t>>>();
    test_MFE<MF<uint16_t, hcd::MontyFullRange<std::uint16_t>>>();
    test_MFE<MF<uint32_t, hcd::MontyFullRange<std::uint32_t>>>();
    test_MFE<MF<uint64_t, hcd::MontyFullRange<std::uint64_t>>>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MFE<MF<__uint128_t, hcd::MontyFullRange<__uint128_t>>>();
#endif
}



TEST(MontgomeryFormExtensions, MontyWrappedStandardMath) {
    namespace hcd = ::hurchalla::detail;
    test_MFE<MF<uint64_t, hcd::MontyWrappedStandardMath<std::uint64_t>>>();
#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
    test_MFE<MF<uint8_t, hcd::MontyWrappedStandardMath<std::uint8_t>>>();
    test_MFE<MF<uint16_t, hcd::MontyWrappedStandardMath<std::uint16_t>>>();
    test_MFE<MF<uint32_t, hcd::MontyWrappedStandardMath<std::uint32_t>>>();
# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MFE<MF<__uint128_t, hcd::MontyWrappedStandardMath<__uint128_t>>>();
# endif
#endif
}

TEST(MontgomeryFormExtensions, MontyFullRangeMasked) {
    namespace hcd = ::hurchalla::detail;
    test_MFE<MF<uint64_t, hcd::MontyFullRangeMasked<std::uint64_t>>>();
#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
    test_MFE<MF<uint8_t, hcd::MontyFullRangeMasked<std::uint8_t>>>();
    test_MFE<MF<uint16_t, hcd::MontyFullRangeMasked<std::uint16_t>>>();
    test_MFE<MF<uint32_t, hcd::MontyFullRangeMasked<std::uint32_t>>>();
# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MFE<MF<__uint128_t,hcd::MontyFullRangeMasked<__uint128_t>>>();
# endif
#endif
}


} // end anonymous namespace
