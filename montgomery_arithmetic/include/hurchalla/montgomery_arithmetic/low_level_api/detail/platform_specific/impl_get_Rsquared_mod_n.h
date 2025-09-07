// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_GET_RSQUARED_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_GET_RSQUARED_MOD_N_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/two_times_restricted.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace detail {


// For discussion purposes, let the unlimited precision constant R represent
// R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t, then
// R = 1<<64.

// Compute (R*R) % n

// Minor note: we use a static member function to disallow ADL.
template <bool nIsGuaranteedLessThanRdiv4, class PTAG>
struct impl_get_Rsquared_mod_n {

  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T n, T inverse_n_modR, T Rmod_n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);

    namespace hc = ::hurchalla;
    T rSquaredModN;
#ifdef HURCHALLA_TESTING_RSQUARED_MOD_N
    if (true) {
#else
    if HURCHALLA_CPP17_CONSTEXPR
            (hc::modular_multiplication_has_slow_perf<T>()) {
#endif
        HPBC_CLOCKWORK_ASSERT2(Rmod_n < n);
        T tmp = Rmod_n;   // Rmod_n ≡ 1*R (mod n)
        int i=0;
        for (; i<8; ++i)
            tmp = hc::modular_addition_prereduced_inputs(tmp, tmp, n);
        // at this point,  tmp ≡ 256*R (mod n)
        constexpr int bitsT = ut_numeric_limits<T>::digits;
        for (; i<bitsT; i*=2) {
            // use montgomery multiplication to square tmp on each iteration
            T u_hi, u_lo;
            u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, tmp, tmp);
            tmp = hc::REDC_standard(u_hi, u_lo, n, inverse_n_modR, PTAG());
        }
        HPBC_CLOCKWORK_ASSERT2(i == bitsT);
        // We should now have  tmp ≡ R*R (mod n).
        // REDC_standard's postcondition guarantees the following:
        HPBC_CLOCKWORK_ASSERT2(tmp < n);

        rSquaredModN = tmp;
        HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN ==
               hc::modular_multiplication_prereduced_inputs(Rmod_n, Rmod_n, n));
    } else {
        rSquaredModN = hc::modular_multiplication_prereduced_inputs(
                                                             Rmod_n, Rmod_n, n);
    }

    HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN < n);
    return rSquaredModN;
  }
};


template<class PTAG>
struct impl_get_Rsquared_mod_n<true, PTAG> {

  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T n, T inverse_n_modR, T Rmod_n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    HPBC_CLOCKWORK_PRECONDITION2(n % 2 == 1);
    HPBC_CLOCKWORK_PRECONDITION2(n > 1);
    // and since the template param nIsGuaranteedLessThanRdiv4 == true,
    constexpr T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 2));
    HPBC_CLOCKWORK_PRECONDITION2(n < Rdiv4);

    namespace hc = ::hurchalla;
    T rSquaredModN;
#ifdef HURCHALLA_TESTING_RSQUARED_MOD_N
    if (true) {
#else
    if HURCHALLA_CPP17_CONSTEXPR
            (hc::modular_multiplication_has_slow_perf<T>()) {
#endif
        HPBC_CLOCKWORK_ASSERT2(Rmod_n < n);
        T tmp = Rmod_n;   // Rmod_n ≡ 1*R (mod n)
        int i=0;

        for (; i<4; ++i)
            tmp = hc::detail::two_times_restricted<T>::call(tmp, n);

        // at this point,  tmp ≡ 16*R (mod n)
        constexpr int bitsT = ut_numeric_limits<T>::digits;

        for (; i<bitsT/2; i*=2) {
            // use montgomery multiplication to square tmp on each iteration
            T u_hi, u_lo;
            u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, tmp, tmp);
            // use the same logic as MontyQuarterRange's montyREDC():
            tmp = hc::REDC_incomplete(u_hi, u_lo, n, inverse_n_modR);
            tmp = static_cast<T>(tmp + n);
            HPBC_CLOCKWORK_ASSERT2(0 < tmp && tmp < static_cast<T>(2*n));
        }
        HPBC_CLOCKWORK_ASSERT2(i == bitsT/2);
        {
            // This final iteration was unrolled from the loop above so we can
            // use standard REDC, which will end with tmp in the range [0, n).
            T u_hi, u_lo;
            u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, tmp, tmp);
            tmp = hc::REDC_standard(u_hi, u_lo, n, inverse_n_modR, PTAG());
        }

        // We should now have  tmp ≡ R*R (mod n).
        // REDC_standard's postcondition guarantees the following:
        HPBC_CLOCKWORK_ASSERT2(tmp < n);

        rSquaredModN = tmp;
        HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN ==
               hc::modular_multiplication_prereduced_inputs(Rmod_n, Rmod_n, n));
    } else {
        rSquaredModN = hc::modular_multiplication_prereduced_inputs(
                                                             Rmod_n, Rmod_n, n);
    }

    HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN < n);
    return rSquaredModN;
  }
};


}} // end namespace

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


#endif
