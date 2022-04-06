// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "test_MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>
#include <cassert>


namespace {


TEST(MontgomeryArithmetic, MontgomeryFormExamples) {
    namespace hc = ::hurchalla;

    // ---- Demonstrate simple modular addition ----
    {
        int64_t modulus = 15;
        int64_t x = 12;
        int64_t y = 4;

        // montgomery arithmetic always requires an odd modulus
        assert(modulus % 2 == 1);
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

        assert(modulus % 2 == 1);  // montgomery arithmetic requires odd modulus
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

TEST(MontgomeryArithmetic, MontyFullRange) {
    namespace hc = ::hurchalla;
    test_custom_monty<hc::detail::MontyFullRange>();
}

TEST(MontgomeryArithmetic, MontyHalfRange) {
    namespace hc = ::hurchalla;
    test_custom_monty<hc::detail::MontyHalfRange>();
}

TEST(MontgomeryArithmetic, MontyQuarterRange) {
    namespace hc = ::hurchalla;
    test_custom_monty<hc::detail::MontyQuarterRange>();
}


TEST(MontgomeryArithmetic, MontyDefault) {
    namespace hc = ::hurchalla;
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

// It would be absolutely normal and expected to use an unsigned integer type
// template argument for MontgomeryForm, but we skip testing them here to save
// on compilation time, because the resulting clases with unsigned int types
// resolve to exactly the same types as will be tested in the MontyFullRange and
// MontyHalfRange and MontyQuarterRange TESTs above.
#if 0
    test_MontgomeryForm<hc::MontgomeryForm<std::uint8_t>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint16_t>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint32_t>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::uint64_t>>();
# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<hc::MontgomeryForm<__uint128_t>>();
# endif
#endif

// To save compilation time we also skip most of the signed integer type tests
// for plain MontgomeryForm.  They should differ from the unsigned versions
// (which in turn map to types we test above) only in the casts they perform for
// convertIn(), convertOut(), max_modulus(), getModulus() and gcd_with_modulus()
// Testing with a single signed integer type should cover much of the potential
// differences regarding casts.
    test_MontgomeryForm<hc::MontgomeryForm<std::int32_t>>();
#if 0
    test_MontgomeryForm<hc::MontgomeryForm<std::int8_t>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::int16_t>>();
    test_MontgomeryForm<hc::MontgomeryForm<std::int64_t>>();
// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_MontgomeryForm<hc::MontgomeryForm<__int128_t>>();
# endif
#endif
}


} // end anonymous namespace
