// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


// We'll define HURCHALLA_ALLOW_INLINE_ASM_ALL here in order to make 
// absolute_value_difference() use an inline asm function version if it is
// available.
// Internally, this inline asm function will also call the generic template
// function version of absolute_value_difference inside a postcondition, in
// order to make sure that the asm result is correct.  Of course postcondition
// checks must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to undefine NDEBUG, which is why we undef
// NDEBUG here too.
#undef HURCHALLA_ALLOW_INLINE_ASM_ALL
#define HURCHALLA_ALLOW_INLINE_ASM_ALL 1
#undef NDEBUG


#include "hurchalla/modular_arithmetic/absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"
#include <cstdint>


template <typename T>
void test_absolute_value_difference()
{
    namespace hc = hurchalla;

    // Test with a few basic examples first
    T a = 5;
    T b = 12;
    EXPECT_TRUE(static_cast<T>(7) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(7) == hc::absolute_value_difference(b, a));
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(b, b));
    a = 7; b = 6;
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(b, a));
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(b, b));

    // --------- Test possible edge cases --------

    a = 0; b = 0;
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(a, b));
    a = 0; b = 1;
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(b, a));
    a = 1; b = 1;
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(a, b));

    a = 0; b = hc::ut_numeric_limits<T>::max();
    EXPECT_TRUE(b == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(b == hc::absolute_value_difference(b, a));
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(b, b));
    a = 1;
    EXPECT_TRUE(static_cast<T>(b-1) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(b-1) == hc::absolute_value_difference(b, a));

    a = 0; b = static_cast<T>(hc::ut_numeric_limits<T>::max() - 1);
    EXPECT_TRUE(b == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(b == hc::absolute_value_difference(b, a));
    a = 1;
    EXPECT_TRUE(static_cast<T>(b-1) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(b-1) == hc::absolute_value_difference(b, a));

    a = static_cast<T>(hc::ut_numeric_limits<T>::max()/2);
    b = static_cast<T>(a + 1);
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(b, a));
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(a, a));
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(b, b));

    b++;
    EXPECT_TRUE(static_cast<T>(2) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(2) == hc::absolute_value_difference(b, a));
    a++;
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(b, a));

    a = static_cast<T>(hc::ut_numeric_limits<T>::max()/2 - 1);
    b = static_cast<T>(a + 1);
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(a, b));
    EXPECT_TRUE(static_cast<T>(1) == hc::absolute_value_difference(b, a));
    EXPECT_TRUE(static_cast<T>(0) == hc::absolute_value_difference(a, a));
}


namespace {
    TEST(ModularArithmetic, absolute_value_difference) {
        test_absolute_value_difference<std::uint8_t>();
        test_absolute_value_difference<std::uint16_t>();
        test_absolute_value_difference<std::uint32_t>();
        test_absolute_value_difference<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_absolute_value_difference<__uint128_t>();
#endif

        test_absolute_value_difference<std::int8_t>();
        test_absolute_value_difference<std::int16_t>();
        test_absolute_value_difference<std::int32_t>();
        test_absolute_value_difference<std::int64_t>();

// It's a slight hack here to use a macro that tells us whether or not the
// compiler supports  __uint128_t, when what we really want is to know is
// whether we can use __int128_t.  Nevertheless in practice, if we have
// __uint128_t then we almost certainly have __int128_t too.
#if HURCHALLA_COMPILER_HAS_UINT128_T()
        test_absolute_value_difference<__int128_t>();
#endif
    }
}
