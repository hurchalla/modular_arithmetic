// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

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
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T u_hi, u_lo;
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rmod_n, a);

    T amodn = static_cast<T>(a % n);
    static_assert(hc::ut_numeric_limits<T>::digits >= 2, "");

    EXPECT_TRUE(hc::REDC_standard(u_hi, u_lo, n, inv_n, hc::LowlatencyTag())
                == amodn);
    EXPECT_TRUE(hc::REDC_standard(u_hi, u_lo, n, inv_n, hc::LowuopsTag())
                == amodn);

    T minuend, subtrahend, result;
    hc::REDC_incomplete(minuend, subtrahend, u_hi, u_lo, n, inv_n, hc::LowlatencyTag());
    result = static_cast<T>(minuend - subtrahend);
    EXPECT_TRUE((minuend < subtrahend) ? static_cast<T>(result+n)==amodn : result==amodn);

    hc::REDC_incomplete(minuend, subtrahend, u_hi, u_lo, n, inv_n, hc::LowuopsTag());
    result = static_cast<T>(minuend - subtrahend);
    EXPECT_TRUE((minuend < subtrahend) ? static_cast<T>(result+n)==amodn : result==amodn);

    T result3 = hc::REDC_incomplete(u_hi, u_lo, n, inv_n, hc::LowlatencyTag());
    EXPECT_TRUE(static_cast<T>(result3+n) == amodn || result3 == amodn);
    T result4 = hc::REDC_incomplete(u_hi, u_lo, n, inv_n, hc::LowuopsTag());
    EXPECT_TRUE(static_cast<T>(result4+n) == amodn || result4 == amodn);
}


template <typename T>
void multi_tests_REDC_identity(T n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

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
template <typename T, class PTAG>
void test_REDCstandard_multiply(T a, T b, T n, T inv_n, T Rsqrd_mod_n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T u_hi, u_lo;
    // convert a and b into montgomery domain
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rsqrd_mod_n, a);
    T a_md = hc::REDC_standard(u_hi, u_lo, n, inv_n, PTAG());
    EXPECT_TRUE(a_md < n);
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rsqrd_mod_n, b);
    T b_md = hc::REDC_standard(u_hi, u_lo, n, inv_n, PTAG());
    EXPECT_TRUE(b_md < n);

    // compute the montgomery domain product of a_md and b_md
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, a_md, b_md);
    T product_md = hc::REDC_standard(u_hi, u_lo, n, inv_n, PTAG());
    EXPECT_TRUE(product_md < n);

    // convert product_md out of montgomery domain, and verify it is correct
    u_hi = 0; u_lo = product_md;
    T product = hc::REDC_standard(u_hi, u_lo, n, inv_n, PTAG());
    T answer = hc::modular_multiplication_prereduced_inputs(
                               static_cast<T>(a % n), static_cast<T>(b % n), n);
    EXPECT_TRUE(product == answer);
}

template <typename T, class PTAG>
void test_REDCincomplete_multiply(T a, T b, T n, T inv_n, T Rsqrd_mod_n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T u_hi, u_lo, minuend, subtrahend;

    // convert a and b into montgomery domain
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rsqrd_mod_n, a);
    hc::REDC_incomplete(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    T a_md = static_cast<T>(minuend - subtrahend);
    a_md = (minuend < subtrahend) ? static_cast<T>(a_md + n) : a_md;
    EXPECT_TRUE(a_md < n);
    T a_md2 = hc::REDC_incomplete(u_hi, u_lo, n, inv_n, PTAG());
    EXPECT_TRUE(a_md == a_md2 || a_md == static_cast<T>(a_md2 + n));

    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, Rsqrd_mod_n, b);
    hc::REDC_incomplete(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    T b_md = static_cast<T>(minuend - subtrahend);
    b_md = (minuend < subtrahend) ? static_cast<T>(b_md + n) : b_md;
    EXPECT_TRUE(b_md < n);
    T b_md2 = hc::REDC_incomplete(u_hi, u_lo, n, inv_n, PTAG());
    EXPECT_TRUE(b_md == b_md2 || b_md == static_cast<T>(b_md2 + n));

    // compute the montgomery domain product of a_md and b_md
    u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, a_md, b_md);
    hc::REDC_incomplete(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    T product_md = static_cast<T>(minuend - subtrahend);
    product_md = (minuend < subtrahend) ? static_cast<T>(product_md + n) : product_md;
    EXPECT_TRUE(product_md < n);
    T product_md2 = hc::REDC_incomplete(u_hi, u_lo, n, inv_n, PTAG());
    EXPECT_TRUE(product_md == product_md2 || product_md == static_cast<T>(product_md2 + n));

    // convert product_md out of montgomery domain, and verify it is correct
    u_hi = 0; u_lo = product_md;
    hc::REDC_incomplete(minuend, subtrahend, u_hi, u_lo, n, inv_n, PTAG());
    T product = static_cast<T>(minuend - subtrahend);
    product = (minuend < subtrahend) ? static_cast<T>(product + n) : product;
    EXPECT_TRUE(product < n);
    T product2 = hc::REDC_incomplete(u_hi, u_lo, n, inv_n, PTAG());
    EXPECT_TRUE(product == product2 || product == static_cast<T>(product2 + n));

    T answer = hc::modular_multiplication_prereduced_inputs(
                               static_cast<T>(a % n), static_cast<T>(b % n), n);
    EXPECT_TRUE(product == answer);
}


template <typename T>
void test_REDC_multiplies(T a, T b, T n, T inv_n, T Rsqrd_mod_n)
{
    test_REDCstandard_multiply<T, hc::LowlatencyTag>(a, b, n, inv_n,
                                                     Rsqrd_mod_n);
    test_REDCstandard_multiply<T, hc::LowuopsTag>(a, b, n, inv_n, Rsqrd_mod_n);
    test_REDCincomplete_multiply<T, hc::LowlatencyTag>(a, b, n, inv_n,
                                                       Rsqrd_mod_n);
    test_REDCincomplete_multiply<T,hc::LowuopsTag>(a, b, n, inv_n, Rsqrd_mod_n);
}


template <typename T>
void multi_tests_REDC_multiply(T n)
{
    static_assert(hc::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(hc::ut_numeric_limits<T>::is_signed), "");
    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    T inv_n = hc::inverse_mod_R(n);
    T Rmod_n = hc::get_R_mod_n(n);
    T R2mod_n = hc::get_Rsquared_mod_n(n, inv_n, Rmod_n);

    // test edge cases, and a couple arbitrary values
    T a = 0; T b = 0;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 0; b = 1;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 1; b = 0;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 1; b = 1;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 2; b = 1;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 1; b = 2;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 2; b = 2;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);

    a = static_cast<T>(static_cast<T>(0) - 1);  b = a;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(static_cast<T>(0) - 1);  b = 1;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(static_cast<T>(0) - 2);
    b = static_cast<T>(static_cast<T>(0) - 1);
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(static_cast<T>(0) - 2);
    b = static_cast<T>(static_cast<T>(0) - 2);
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);

    a = n; b = 5;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = n; b = n;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(n - 1); b = 3;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(n - 1); b = a;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>(n + 1);   // this might overflow to 0, which is ok if so.
    b = static_cast<T>(n - 1);
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);

    a = static_cast<T>(n/2); b = a;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = static_cast<T>((n/2)+1); b = static_cast<T>(n/2);
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);

    a = 127; b = 13;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 200; b = 254;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
    a = 93; b = 12;
    test_REDC_multiplies(a, b, n, inv_n, R2mod_n);
}



template <typename T>
void REDC_test_all(T n)
{
    multi_tests_REDC_identity(n);
    multi_tests_REDC_multiply(n);
}



} // end unnamed namespace

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
