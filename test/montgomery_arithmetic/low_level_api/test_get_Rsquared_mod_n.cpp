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


using namespace ::hurchalla;


template <typename T>
void test_single_R2(T n)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    constexpr int digitsR = ut_numeric_limits<T>::digits;
    constexpr T Rdiv4 = static_cast<T>(1) << (digitsR - 2);

    T rmodn = get_R_mod_n(n);
    T inv = inverse_mod_R(n);
    using P = typename safely_promote_unsigned<T>::type;
    T one = static_cast<T>(1);
    // the next line tests inverse_mod_R - we might as well test it while here.
    EXPECT_TRUE(static_cast<T>(static_cast<P>(inv) * static_cast<P>(n)) == one);
    T answer = modular_multiplication_prereduced_inputs(rmodn, rmodn, n);

    if (n < Rdiv4) {
        T r2modn_1 = get_Rsquared_mod_n<T, true, LowlatencyTag>(n, inv, rmodn);
        T r2modn_2 = get_Rsquared_mod_n<T, true, LowuopsTag>(n, inv, rmodn);
        EXPECT_TRUE(r2modn_1 == answer);
        EXPECT_TRUE(r2modn_2 == answer);
    }
    // test version that works for all n
    {
        T r2modn_1 = get_Rsquared_mod_n<T, false, LowlatencyTag>(n, inv, rmodn);
        T r2modn_2 = get_Rsquared_mod_n<T, false, LowuopsTag>(n, inv, rmodn);
        EXPECT_TRUE(r2modn_1 == answer);
        EXPECT_TRUE(r2modn_2 == answer);
    }
}


template <std::size_t ARRAY_SIZE, typename T>
void test_single_R2_array(T n)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    constexpr int digitsR = ut_numeric_limits<T>::digits;
    constexpr T Rdiv4 = static_cast<T>(1) << (digitsR - 2);

    std::array<T, ARRAY_SIZE> a_n;
    std::array<T, ARRAY_SIZE> a_rmn;
    std::array<T, ARRAY_SIZE> a_invn;
    std::array<T, ARRAY_SIZE> answer;
    for (std::size_t i=0; i<ARRAY_SIZE; ++i) {
        if (n < 3 + 2*i)
            a_n[i] = 3;
        else
            a_n[i] = static_cast<T>(n - 2*i);
        a_rmn[i] = get_R_mod_n(a_n[i]);
        a_invn[i] = inverse_mod_R(a_n[i]);
        answer[i] = modular_multiplication_prereduced_inputs(
                                                    a_rmn[i], a_rmn[i], a_n[i]);
        // we might as well test inverse_mod_R while here.
        using P = typename safely_promote_unsigned<T>::type;
        T one = static_cast<T>(1);
        EXPECT_TRUE(static_cast<T>(
                    static_cast<P>(a_invn[i]) * static_cast<P>(a_n[i])) == one);
    }

    // since we subtracted from n to set a_n, (a_n[0] < Rdiv4) covers all a_n[i]
    if (a_n[0] < Rdiv4) {
        auto r2mn1 = get_Rsquared_mod_n<T, ARRAY_SIZE, true, LowlatencyTag>
                                                           (a_n, a_invn, a_rmn);
        auto r2mn2 = get_Rsquared_mod_n<T, ARRAY_SIZE, true, LowuopsTag>
                                                           (a_n, a_invn, a_rmn);
        EXPECT_TRUE(r2mn1 == answer);
        EXPECT_TRUE(r2mn2 == answer);
    }
    // test version that works for any size a_n[i]
    {
        auto r2mn1 = get_Rsquared_mod_n<T, ARRAY_SIZE, false, LowlatencyTag>
                                                           (a_n, a_invn, a_rmn);
        auto r2mn2 = get_Rsquared_mod_n<T, ARRAY_SIZE, false, LowuopsTag>
                                                           (a_n, a_invn, a_rmn);
        EXPECT_TRUE(r2mn1 == answer);
        EXPECT_TRUE(r2mn2 == answer);
    }
}



template <typename T>
void test_R2_exhaustive()
{
    T max = ut_numeric_limits<T>::max();
    EXPECT_TRUE(max > 0);
    T evenmax = static_cast<T>((max/2)*2);
    T oddmax = (evenmax != max) ? max : static_cast<T>(max - 1);
    // get_Rsquared_mod_n's preconditions require input n is odd and > 1.
    for (T n=oddmax; n>1; n=static_cast<T>(n-2)) {
        test_single_R2(n);
        test_single_R2_array<3>(n);  // array size of 3 is an arbitrary size
    }
}


template <typename T>
void test_R2()
{
    T max = ut_numeric_limits<T>::max();
    EXPECT_TRUE(max > 0);
    T evenmax = static_cast<T>((max/2)*2);
    T oddmax = (evenmax != max) ? max : static_cast<T>(max - 1);
    T oddquartermax = static_cast<T>((max/8)*2 + 1);

    // get_Rsquared_mod_n's preconditions require input n is odd and > 1.

    T n = 3;
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = 9;
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = 11;
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = 21;
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = oddmax;
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = static_cast<T>(oddmax - 2);
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = static_cast<T>(oddmax - 6);
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = oddquartermax;
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = static_cast<T>(oddquartermax + 2);
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);

    n = static_cast<T>(oddquartermax - 2);
    test_single_R2(n);
    test_single_R2_array<1>(n);
    test_single_R2_array<2>(n);
    test_single_R2_array<5>(n);
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
