// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace {


namespace hc = ::hurchalla;

// verify that  REDC(a*(R mod n)) == a
template <typename T>
void test_REDC_identity(T a, T n, T inv_n, T Rmod_n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T u_hi, u_lo;
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rmod_n, a);

    T amodn = static_cast<T>(a % n);
    static_assert(hc::ut_numeric_limits<T>::digits >= 2, "");
    T Rdiv4 = static_cast<T>(
                   static_cast<T>(1) << (hc::ut_numeric_limits<T>::digits - 2));
    {  // Fullrange has no requirements on n
        EXPECT_TRUE(hc::REDC(u_hi, u_lo, n, inv_n) == amodn);
        EXPECT_TRUE(hc::REDC(u_hi, u_lo, n, inv_n, hc::FullrangeTag(),
                               hc::LowlatencyTag()) == amodn);
        EXPECT_TRUE(hc::REDC(u_hi, u_lo, n, inv_n, hc::FullrangeTag(),
                               hc::LowuopsTag()) == amodn);
    }
    if (n < Rdiv4) {
        T result1 = hc::REDC(u_hi, u_lo, n, inv_n, hc::QuarterrangeTag(),
                               hc::LowlatencyTag());
        T result2 = hc::REDC(u_hi, u_lo, n, inv_n, hc::QuarterrangeTag(),
                               hc::LowuopsTag());
        EXPECT_TRUE(result1 == amodn || result1 == amodn + n);
        EXPECT_TRUE(result2 == amodn || result2 == amodn + n);
    }
}


template <typename T>
void multi_tests_REDC_identity(T n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T inv_n = hc::inverse_mod_R(n);
    T Rmod_n = hc::get_R_mod_n(n);

    // test edge cases, and a couple arbitrary values
    T a = 0;
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = 1;
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = 2;
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = static_cast<T>(static_cast<T>(0) - 1);
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = static_cast<T>(static_cast<T>(0) - 2);
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = n;
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = static_cast<T>(n - 1);
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = static_cast<T>(n + 1);   // this might overflow to 0, which is ok if so.
    test_REDC_identity(a, n, inv_n, Rmod_n);

    a = static_cast<T>(n/2);
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = static_cast<T>((n/2)+1);
    test_REDC_identity(a, n, inv_n, Rmod_n);

    a = 127;
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = 200;
    test_REDC_identity(a, n, inv_n, Rmod_n);
    a = 93;
    test_REDC_identity(a, n, inv_n, Rmod_n);
}


// verify that montgomery multiplication (via REDC) provides the same answer
// as standard modular multiplication
template <typename T, class MTAG, class PTAG>
void test_REDC_multiply(T a, T b, T n, T inv_n, T Rsqrd_mod_n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    T u_hi, u_lo;
    // convert a and b into montgomery domain
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rsqrd_mod_n, a);
    T a_md = hc::REDC(u_hi, u_lo, n, inv_n, MTAG(), PTAG());
    EXPECT_TRUE((std::is_same<MTAG, hc::FullrangeTag>::value) ?
                 a_md < n : a_md < 2*n);
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rsqrd_mod_n, b);
    T b_md = hc::REDC(u_hi, u_lo, n, inv_n, MTAG(), PTAG());
    EXPECT_TRUE((std::is_same<MTAG, hc::FullrangeTag>::value) ?
                 b_md < n : b_md < 2*n);

    // compute the montgomery domain product of a_md and b_md
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, a_md, b_md);
    T product_md = hc::REDC(u_hi, u_lo, n, inv_n, MTAG(), PTAG());
    EXPECT_TRUE((std::is_same<MTAG, hc::FullrangeTag>::value) ?
                 product_md < n : product_md < 2*n);

    // convert product_md out of montgomery domain, and verify it is correct
    u_hi = 0; u_lo = product_md;
    T product = hc::REDC(u_hi, u_lo, n, inv_n, MTAG(), PTAG());
    T answer = hc::modular_multiplication_prereduced_inputs(
                               static_cast<T>(a % n), static_cast<T>(b % n), n);
    EXPECT_TRUE((std::is_same<MTAG, hc::FullrangeTag>::value) ?
                 product == answer :
                 product == answer || product == answer + n);
}


template <typename T, class MTAG, class PTAG>
void multi_tests_REDC_multiply(T n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    static_assert(hc::ut_numeric_limits<T>::digits >= 2, "");
    T Rdiv4 = static_cast<T>(
                   static_cast<T>(1) << (hc::ut_numeric_limits<T>::digits - 2));

    // if n is invalid for the MTAG, change n to be some value we can test
    if (std::is_same<MTAG, hc::QuarterrangeTag>::value) {
        while (n >= Rdiv4)
            n = static_cast<T>(n - Rdiv4);
    }
    if (n < 3)
        n = 3;
    if (n%2 != 1)
        --n;

    T inv_n = hc::inverse_mod_R(n);
    T Rmod_n = hc::get_R_mod_n(n);
    T R2mod_n = hc::get_Rsquared_mod_n(n, inv_n, Rmod_n);

    // test edge cases, and a couple arbitrary values
    T a = 0; T b = 0;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 0; b = 1;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 1; b = 0;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 1; b = 1;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 2; b = 1;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 1; b = 2;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 2; b = 2;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);

    a = static_cast<T>(static_cast<T>(0) - 1);  b = a;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(static_cast<T>(0) - 1);  b = 1;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(static_cast<T>(0) - 2);
    b = static_cast<T>(static_cast<T>(0) - 1);
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(static_cast<T>(0) - 2);
    b = static_cast<T>(static_cast<T>(0) - 2);
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);

    a = n; b = 5;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = n; b = n;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(n - 1); b = 3;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(n - 1); b = a;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(n + 1);   // this might overflow to 0, which is ok if so.
    b = static_cast<T>(n - 1);
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);

    a = static_cast<T>(n/2); b = a;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>((n/2)+1); b = static_cast<T>(n/2);
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);

    a = 127; b = 13;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 200; b = 254;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
    a = 93; b = 12;
    test_REDC_multiply<T, MTAG, PTAG>(a, b, n, inv_n, R2mod_n);
}


template <typename T, class MTAG, class PTAG>
void test_REDC_is_zero(T a, T n, T inv_n, T Rmod_n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    // get the reference answer on whether the result should be zero.
    T amodn = static_cast<T>(a % n);
    bool isZeroAnswer = (amodn == 0);

    // test that the REDC resultIsZero matches what we expect
    T u_hi, u_lo;
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rmod_n, a);
    bool resultIsZero;
    T result = REDC(u_hi, u_lo, n, inv_n, resultIsZero, MTAG(), PTAG());
    EXPECT_TRUE(resultIsZero == isZeroAnswer);

    // we might as well test that the result itself is correct too,
    // even though it's not our focus here.
    EXPECT_TRUE((std::is_same<MTAG, hc::FullrangeTag>::value) ?
                 result == amodn :
                 result == amodn || result == amodn + n);
}


template <typename T, class MTAG, class PTAG>
void multi_tests_REDC_is_zero(T n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    static_assert(hc::ut_numeric_limits<T>::digits >= 2, "");
    T Rdiv4 = static_cast<T>(
                   static_cast<T>(1) << (hc::ut_numeric_limits<T>::digits - 2));

    // if n is invalid for the MTAG, change n to be some value we can test
    if (std::is_same<MTAG, hc::QuarterrangeTag>::value) {
        while (n >= Rdiv4)
            n = static_cast<T>(n - Rdiv4);
    }
    if (n < 3)
        n = 3;
    if (n%2 != 1)
        --n;

    T inv_n = hc::inverse_mod_R(n);
    T Rmod_n = hc::get_R_mod_n(n);

    T a = 0;
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    a = 1;
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    a = n;
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    a = static_cast<T>(2*n);  // this might overflow, but that's ok
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    ++a;
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    a = static_cast<T>(a - 2);
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    a = static_cast<T>(static_cast<T>(0) - Rmod_n);
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    ++a;
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
    a = static_cast<T>(a - 2);
    test_REDC_is_zero<T, MTAG, PTAG>(a, n, inv_n, Rmod_n);
}


template <typename T>
void REDC_test_all(T n)
{
    multi_tests_REDC_identity(n);

    multi_tests_REDC_multiply<T, hc::FullrangeTag, hc::LowlatencyTag>(n);
    multi_tests_REDC_multiply<T, hc::FullrangeTag, hc::LowuopsTag>(n);
    multi_tests_REDC_multiply<T, hc::QuarterrangeTag, hc::LowlatencyTag>(n);
    multi_tests_REDC_multiply<T, hc::QuarterrangeTag, hc::LowuopsTag>(n);

    multi_tests_REDC_is_zero<T, hc::FullrangeTag, hc::LowlatencyTag>(n);
    multi_tests_REDC_is_zero<T, hc::FullrangeTag, hc::LowuopsTag>(n);
    multi_tests_REDC_is_zero<T, hc::QuarterrangeTag, hc::LowlatencyTag>(n);
    multi_tests_REDC_is_zero<T, hc::QuarterrangeTag, hc::LowuopsTag>(n);
}



TEST(MontgomeryArithmetic, REDC8) {
    std::vector<uint8_t> moduli { 3, 255, 19, 21, 211, 23, 171 };
    for (auto n : moduli)
        REDC_test_all(n);
}
TEST(MontgomeryArithmetic, REDC16) {
    std::vector<uint16_t> moduli { 3, 17, UINT16_C(65535),
                          UINT16_C(65533), UINT16_C(357), UINT16_C(32253),
                          UINT16_C(11111) };
    for (auto n : moduli)
        REDC_test_all(n);
}
TEST(MontgomeryArithmetic, REDC32) {
    std::vector<uint32_t> moduli { 3, 13, UINT32_C(4294967295),
                          UINT32_C(4294967293), UINT32_C(2147483347),
                          UINT32_C(246098243), UINT32_C(1111111) };
    for (auto n : moduli)
        REDC_test_all(n);
}
TEST(MontgomeryArithmetic, REDC64) {
    std::vector<uint64_t> moduli { 3, 11, UINT64_C(18446744073709551615),
                          UINT64_C(18446744073709551613),
                          UINT64_C(4294967295),
                          UINT64_C(3194806714689), UINT64_C(11111111311) };
    for (auto n : moduli)
        REDC_test_all(n);
}

#if !defined(__GNUC__) || __GNUC__ >= 11 || defined(__INTEL_COMPILER) || \
                                            defined(__clang__)
// GCC has a compiler bug that causes an incorrect value of n to be produced and
// thus results in one of my google test assertions failing.  See
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98474 .  The bug appears to have
// been introduced as a regression to gcc in v5.1.  It exists up to the latest
// released version (v10.2) of gcc at the time of this writing.  It's unclear
// at the moment whether __uint128_t is safe to use with any version of gcc
// between 5.1 and 10.2.  The patch appears to fix the bug, and it is scheduled
// to be in the gcc 11 release.
// For now, the #if above disables the following tests on gcc, since they will
// fail on gcc at optimization level -O1 or higher due to the compiler bug.
# if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(MontgomeryArithmetic, REDC128) {
    __uint128_t zero = 0;
    std::vector<__uint128_t> moduli { 3, 11, zero-1, zero-3,
                  static_cast<__uint128_t>(UINT64_C(18446744073709551613)) *
                                             UINT64_C(18446744073709551611),
                  static_cast<__uint128_t>(UINT64_C(35698723439051265)) *
                                                UINT64_C(70945870135873583),
                  static_cast<__uint128_t>(UINT64_C(34069834503)) *
                                              UINT64_C(895835939) };
    for (auto n : moduli)
        REDC_test_all(n);
}
# endif
#endif


} // end unnamed namespace

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
