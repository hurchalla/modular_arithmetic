// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

// lets us get inside the black box of get_Rsquared_mod_n() to ensure that
// we test the complex compiled possibility rather than the trivial one.
#define HURCHALLA_TESTING_RSQUARED_MOD_N 1

#include "hurchalla/montgomery_arithmetic/low_level_api/get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>

namespace {


namespace hc = ::hurchalla;


template <typename T>
void test_single_R2(T n)
{
    using P = typename hc::safely_promote_unsigned<T>::type;
    T one = static_cast<T>(1);
    T inv = hc::inverse_mod_R(n);
    // the next line tests inverse_mod_R - we might as well test it while here.
    EXPECT_TRUE(static_cast<T>(static_cast<P>(inv) * static_cast<P>(n)) == one);

    T rmodn = hc::get_R_mod_n(n);
    T r2modn_1 = hc::get_Rsquared_mod_n(n, inv, rmodn);
    T r2modn_2 = hc::modular_multiplication_prereduced_inputs(rmodn, rmodn, n);
    // finally, test get_Rsquared_mod_n()
    EXPECT_TRUE(r2modn_1 == r2modn_2);
}


template <typename T>
void test_R2_exhaustive()
{
    T max = hc::ut_numeric_limits<T>::max();
    EXPECT_TRUE(max > 0);
    T evenmax = static_cast<T>((max/2)*2);
    T oddmax = (evenmax != max) ? max : static_cast<T>(max - 1);
    // get_Rsquared_mod_n's preconditions require input n is odd and > 1.
    for (T n=oddmax; n>1; n=static_cast<T>(n-2))
        test_single_R2(n);
}


template <typename T>
void test_R2()
{
    T max = hc::ut_numeric_limits<T>::max();
    EXPECT_TRUE(max > 0);
    T evenmax = static_cast<T>((max/2)*2);
    T oddmax = (evenmax != max) ? max : static_cast<T>(max - 1);
    T oddhalfmax = static_cast<T>((max/4)*2 + 1);

    // get_Rsquared_mod_n's preconditions require input n is odd and > 1.

    test_single_R2(static_cast<T>(3));
    test_single_R2(static_cast<T>(9));
    test_single_R2(static_cast<T>(11));
    test_single_R2(static_cast<T>(21));

    test_single_R2(static_cast<T>(oddmax));
    test_single_R2(static_cast<T>(oddmax - 2));
    test_single_R2(static_cast<T>(oddmax - 6));

    test_single_R2(static_cast<T>(oddhalfmax));
    test_single_R2(static_cast<T>(oddhalfmax + 2));
    test_single_R2(static_cast<T>(oddhalfmax - 2));
}



TEST(MontgomeryArithmetic, get_Rsquared_mod_N) {
    test_R2<std::uint8_t>();
    test_R2<std::uint16_t>();
    test_R2<std::uint32_t>();
    test_R2<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_R2<__uint128_t>();
#endif

    test_R2_exhaustive<std::uint8_t>();
    test_R2_exhaustive<std::uint16_t>();
}


} // end unnamed namespace
