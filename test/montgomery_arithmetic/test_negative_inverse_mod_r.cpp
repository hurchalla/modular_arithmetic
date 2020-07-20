// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


#include "hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>


template <typename T>
void test_single_inverse(T a)
{
    namespace mont = hurchalla::montgomery_arithmetic;
    using U = typename mont::safely_promote_unsigned<T>::type;
    T minusOne = static_cast<T>(static_cast<U>(0) - static_cast<U>(1));

    T inv = mont::negative_inverse_mod_r(a);
    EXPECT_TRUE(static_cast<T>(static_cast<U>(inv) * static_cast<U>(a)) ==
                                                                      minusOne);
}


template <typename T>
void test_inverse_exhaustive()
{
    namespace ma = hurchalla::modular_arithmetic;
    T tmax = ma::ma_numeric_limits<T>::max();
    T evenmax = static_cast<T>((tmax/2)*2);
    T oddmax = (evenmax != tmax) ? tmax : static_cast<T>(tmax - 1);

    for (T a=oddmax; a>1; a=static_cast<T>(a-2))
        test_single_inverse(a);
    test_single_inverse(static_cast<T>(1));
}


template <typename T>
void test_negative_inverse_mod_r()
{
    namespace ma = hurchalla::modular_arithmetic;
    T tmax = ma::ma_numeric_limits<T>::max();
    T evenmax = static_cast<T>((tmax/2)*2);
    T oddmax = (evenmax != tmax) ? tmax : static_cast<T>(tmax - 1);
    T oddhalfmax = static_cast<T>((tmax/4)*2 + 1);

    // negative_inverse_mod_r's preconditions require input a is odd.

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
    TEST(MontgomeryArithmetic, negative_inverse_mod_r) {
        test_negative_inverse_mod_r<uint8_t>();
        test_negative_inverse_mod_r<uint16_t>();
        test_negative_inverse_mod_r<uint32_t>();
        test_negative_inverse_mod_r<uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_negative_inverse_mod_r<__uint128_t>();
#endif

        test_inverse_exhaustive<uint8_t>();
        test_inverse_exhaustive<uint16_t>();
    }
}
