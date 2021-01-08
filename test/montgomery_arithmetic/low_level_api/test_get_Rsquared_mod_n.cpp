// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

// lets us get inside the black box of get_Rsquared_mod_n() to ensure that
// we test the complex compiled possibility rather than the trivial one.
#define HURCHALLA_TESTING_RSQUARED_MOD_N 1

#include "hurchalla/montgomery_arithmetic/low_level_api/get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>


template <typename T, class MTAG>
T get_max_allowable_modulus()
{
    namespace mont = hurchalla::montgomery_arithmetic;
    namespace ut = hurchalla::util;
    T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                        (ut::ut_numeric_limits<T>::digits - 1));
    T Rdiv4 = static_cast<T>(Rdiv2 / 2);
    T Rdiv6 = static_cast<T>(Rdiv2 / 3);

    if (std::is_same<MTAG, mont::FullrangeTag>::value)
        return static_cast<T>(Rdiv2 - 1 + Rdiv2);
    else if (std::is_same<MTAG, mont::HalfrangeTag>::value)
        return static_cast<T>(Rdiv2 - 1);
    else if (std::is_same<MTAG, mont::QuarterrangeTag>::value)
        return static_cast<T>(Rdiv4 - 1);
    else if (std::is_same<MTAG, mont::SixthrangeTag>::value)
        return static_cast<T>(Rdiv6 - 1);
    else {
        EXPECT_TRUE(false);  // there should be no other tag types than above
        return static_cast<T>(0);
    }
}


template <typename T, class MTAG>
void test_single_R2(T n)
{
    namespace ma = hurchalla::modular_arithmetic;
    namespace mont = hurchalla::montgomery_arithmetic;
    namespace ut = hurchalla::util;
    using P = typename ut::safely_promote_unsigned<T>::type;
    // a failure on the next line would mean the test case was written wrong.
    T max = get_max_allowable_modulus<T, MTAG>();
    EXPECT_TRUE(n <= max);
    T one = static_cast<T>(1);

    T inv = mont::inverse_mod_R(n);
    // the next line tests inverse_mod_R - we might as well test it while here.
    EXPECT_TRUE(static_cast<T>(static_cast<P>(inv) * static_cast<P>(n)) == one);

    T rmodn = mont::get_R_mod_n(n);
    T r2modn_1 = mont::get_Rsquared_mod_n(n, inv, rmodn, MTAG());
    T r2modn_2 = ma::modular_multiplication_prereduced_inputs(rmodn, rmodn, n);
    // finally, test get_Rsquared_mod_n()
    EXPECT_TRUE(r2modn_1 == r2modn_2);
}


template <typename T, class MTAG>
void test_R2_exhaustive()
{
    T max = get_max_allowable_modulus<T, MTAG>();
    T evenmax = static_cast<T>((max/2)*2);
    T oddmax = (evenmax != max) ? max : static_cast<T>(max - 1);

    // get_Rsquared_mod_n's preconditions require input n is odd and > 1.
    for (T n=oddmax; n>1; n=static_cast<T>(n-2))
        test_single_R2<T, MTAG>(n);
}


template <typename T, class MTAG>
void test_R2()
{
    T max = get_max_allowable_modulus<T, MTAG>();
    T evenmax = static_cast<T>((max/2)*2);
    T oddmax = (evenmax != max) ? max : static_cast<T>(max - 1);
    T oddhalfmax = static_cast<T>((max/4)*2 + 1);

    // get_Rsquared_mod_n's preconditions require input n is odd and > 1.

    test_single_R2<T, MTAG>(static_cast<T>(3));
    test_single_R2<T, MTAG>(static_cast<T>(9));
    test_single_R2<T, MTAG>(static_cast<T>(11));
    test_single_R2<T, MTAG>(static_cast<T>(21));

    test_single_R2<T, MTAG>(static_cast<T>(oddmax));
    test_single_R2<T, MTAG>(static_cast<T>(oddmax - 2));
    test_single_R2<T, MTAG>(static_cast<T>(oddmax - 6));

    test_single_R2<T, MTAG>(static_cast<T>(oddhalfmax));
    test_single_R2<T, MTAG>(static_cast<T>(oddhalfmax + 2));
    test_single_R2<T, MTAG>(static_cast<T>(oddhalfmax - 2));
}


namespace {
    TEST(MontgomeryArithmetic, get_Rsquared_mod_N) {
        namespace mont = hurchalla::montgomery_arithmetic;

        test_R2<std::uint8_t, mont::FullrangeTag>();
        test_R2<std::uint8_t, mont::HalfrangeTag>();
        test_R2<std::uint8_t, mont::QuarterrangeTag>();
        test_R2<std::uint8_t, mont::SixthrangeTag>();

        test_R2<std::uint16_t, mont::FullrangeTag>();
        test_R2<std::uint16_t, mont::HalfrangeTag>();
        test_R2<std::uint16_t, mont::QuarterrangeTag>();
        test_R2<std::uint16_t, mont::SixthrangeTag>();

        test_R2<std::uint32_t, mont::FullrangeTag>();
        test_R2<std::uint32_t, mont::HalfrangeTag>();
        test_R2<std::uint32_t, mont::QuarterrangeTag>();
        test_R2<std::uint32_t, mont::SixthrangeTag>();

        test_R2<std::uint64_t, mont::FullrangeTag>();
        test_R2<std::uint64_t, mont::HalfrangeTag>();
        test_R2<std::uint64_t, mont::QuarterrangeTag>();
        test_R2<std::uint64_t, mont::SixthrangeTag>();

#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_R2<__uint128_t, mont::FullrangeTag>();
        test_R2<__uint128_t, mont::HalfrangeTag>();
        test_R2<__uint128_t, mont::QuarterrangeTag>();
        test_R2<__uint128_t, mont::SixthrangeTag>();
#endif

        test_R2_exhaustive<std::uint8_t, mont::FullrangeTag>();
        test_R2_exhaustive<std::uint8_t, mont::HalfrangeTag>();
        test_R2_exhaustive<std::uint8_t, mont::QuarterrangeTag>();
        test_R2_exhaustive<std::uint8_t, mont::SixthrangeTag>();

        test_R2_exhaustive<std::uint16_t, mont::FullrangeTag>();
        test_R2_exhaustive<std::uint16_t, mont::HalfrangeTag>();
        test_R2_exhaustive<std::uint16_t, mont::QuarterrangeTag>();
        test_R2_exhaustive<std::uint16_t, mont::SixthrangeTag>();
    }
}
