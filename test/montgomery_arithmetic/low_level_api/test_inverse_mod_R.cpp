// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>


template <typename T>
void test_single_inverse(T a)
{
    namespace hc = hurchalla;
    using P = typename hc::safely_promote_unsigned<T>::type;
    T one = static_cast<T>(1);

    T inv = hc::inverse_mod_R(a);
    EXPECT_TRUE(static_cast<T>(static_cast<P>(inv) * static_cast<P>(a)) == one);
    // the #if is a slight hack, but inverse_mod_R is only constexpr for C++14
    // and above (C++11's support for constexpr functions was too primitive)
#if (__cplusplus >= 201402L) || \
        (defined(_MSVC_LANG) && _MSVC_LANG >= 201402L && _MSC_VER >= 1910)
    // test constexpr use of inverse_mod_R
    static_assert(hc::inverse_mod_R(static_cast<T>(1)) == 1, "");
#endif
}


template <typename T>
void test_inverse_exhaustive()
{
    namespace hc = hurchalla;
    T tmax = hc::ut_numeric_limits<T>::max();
    T evenmax = static_cast<T>((tmax/2)*2);
    T oddmax = (evenmax != tmax) ? tmax : static_cast<T>(tmax - 1);

    for (T a=oddmax; a>1; a=static_cast<T>(a-2))
        test_single_inverse(a);
    test_single_inverse(static_cast<T>(1));
}


template <typename T>
void test_inverse_mod_r()
{
    namespace hc = hurchalla;
    T tmax = hc::ut_numeric_limits<T>::max();
    T evenmax = static_cast<T>((tmax/2)*2);
    T oddmax = (evenmax != tmax) ? tmax : static_cast<T>(tmax - 1);
    T oddhalfmax = static_cast<T>((tmax/4)*2 + 1);

    // inverse_mod_r's preconditions require input a is odd.

    test_single_inverse(static_cast<T>(1));
    test_single_inverse(static_cast<T>(3));
    test_single_inverse(static_cast<T>(5));
    test_single_inverse(static_cast<T>(7));

    test_single_inverse(static_cast<T>(oddmax));
    test_single_inverse(static_cast<T>(oddmax - 2));
    test_single_inverse(static_cast<T>(oddmax - 4));

    test_single_inverse(static_cast<T>(oddhalfmax));
    test_single_inverse(static_cast<T>(oddhalfmax + 2));
    test_single_inverse(static_cast<T>(oddhalfmax - 2));
}


namespace {
    TEST(MontgomeryArithmetic, inverse_mod_r) {
        test_inverse_mod_r<std::uint8_t>();
        test_inverse_mod_r<std::uint16_t>();
        test_inverse_mod_r<std::uint32_t>();
        test_inverse_mod_r<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_inverse_mod_r<__uint128_t>();
#endif

        test_inverse_exhaustive<std::uint8_t>();
        test_inverse_exhaustive<std::uint16_t>();
    }
}
