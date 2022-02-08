// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


#include "hurchalla/montgomery_arithmetic/detail/experimental/negative_inverse_mod_R.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>

namespace {


namespace hc = ::hurchalla;

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4309)
#endif

template <typename T>
void test_single_negative_inverse(T a)
{
    using P = typename hc::safely_promote_unsigned<T>::type;
    T minusOne = static_cast<T>(static_cast<P>(0) - static_cast<P>(1));

    T inv = hc::detail::negative_inverse_mod_R(a);
    EXPECT_TRUE(static_cast<T>(static_cast<P>(inv) * static_cast<P>(a)) ==
                                                                      minusOne);
}

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


template <typename T>
void test_negative_inverse_exhaustive()
{
    T tmax = hc::ut_numeric_limits<T>::max();
    T evenmax = static_cast<T>((tmax/2)*2);
    T oddmax = (evenmax != tmax) ? tmax : static_cast<T>(tmax - 1);

    for (T a=oddmax; a>1; a=static_cast<T>(a-2))
        test_single_negative_inverse(a);
    test_single_negative_inverse(static_cast<T>(1));
}


template <typename T>
void test_negative_inverse_mod_R()
{
    T tmax = hc::ut_numeric_limits<T>::max();
    T evenmax = static_cast<T>((tmax/2)*2);
    T oddmax = (evenmax != tmax) ? tmax : static_cast<T>(tmax - 1);
    T oddhalfmax = static_cast<T>((tmax/4)*2 + 1);

    // negative_inverse_mod_r's preconditions require input a is odd.

    test_single_negative_inverse(static_cast<T>(1));
    test_single_negative_inverse(static_cast<T>(3));
    test_single_negative_inverse(static_cast<T>(5));
    test_single_negative_inverse(static_cast<T>(7));

    test_single_negative_inverse(static_cast<T>(oddmax));
    test_single_negative_inverse(static_cast<T>(oddmax - 2));
    test_single_negative_inverse(static_cast<T>(oddmax - 4));

    test_single_negative_inverse(static_cast<T>(oddhalfmax));
    test_single_negative_inverse(static_cast<T>(oddhalfmax + 2));
    test_single_negative_inverse(static_cast<T>(oddhalfmax - 2));
}



TEST(MontgomeryArithmetic, negative_inverse_mod_r) {
    test_negative_inverse_mod_R<std::uint8_t>();
    test_negative_inverse_mod_R<std::uint16_t>();
    test_negative_inverse_mod_R<std::uint32_t>();
    test_negative_inverse_mod_R<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_negative_inverse_mod_R<__uint128_t>();
#endif

    test_negative_inverse_exhaustive<std::uint8_t>();
    test_negative_inverse_exhaustive<std::uint16_t>();
}


} // end unnamed namespace
