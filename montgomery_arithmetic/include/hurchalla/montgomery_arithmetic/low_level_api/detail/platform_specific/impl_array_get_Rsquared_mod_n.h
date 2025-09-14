// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_ARRAY_GET_RSQUARED_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_ARRAY_GET_RSQUARED_MOD_N_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/two_times_restricted.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <array>

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
struct impl_array_get_Rsquared_mod_n {

  template <typename T, std::size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE static
  std::array<T, ARRAY_SIZE> call(const std::array<T, ARRAY_SIZE>& n,
                                 const std::array<T, ARRAY_SIZE>& inverse_n_modR,
                                 const std::array<T, ARRAY_SIZE>& Rmod_n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    using std::size_t;
    namespace hc = ::hurchalla;
    if (HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE) {
        for (size_t j=0; j<ARRAY_SIZE; ++j) {
            HPBC_CLOCKWORK_PRECONDITION2(n[j] % 2 == 1);
            HPBC_CLOCKWORK_PRECONDITION2(n[j] > 1);
            HPBC_CLOCKWORK_PRECONDITION2(Rmod_n[j] < n[j]);
        }
    }
    std::array<T, ARRAY_SIZE> rSquaredModN;

#ifdef HURCHALLA_TESTING_RSQUARED_MOD_N
    if (true) {
#else
    if (hc::modular_multiplication_has_slow_perf<T>()) {
#endif
        int i=0;
        std::array<T, ARRAY_SIZE> tmp(Rmod_n);   // Rmod_n ≡ 1*R (mod n)
        for (; i<8; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                tmp[j] = hc::modular_addition_prereduced_inputs(tmp[j], tmp[j], n[j]);
        }
        // at this point,  tmp ≡ 256*R (mod n)
        constexpr int bitsT = ut_numeric_limits<T>::digits;

        for (; i<bitsT; i*=2) {
            // use montgomery multiplication to square tmp on each iteration
            T u_hi, u_lo;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, tmp[j], tmp[j]);
                tmp[j] = hc::REDC_standard(u_hi, u_lo, n[j], inverse_n_modR[j], PTAG());
            }
        }
        HPBC_CLOCKWORK_ASSERT2(i == bitsT);
        rSquaredModN = tmp;

        if (HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE) {
            for (size_t j=0; j<ARRAY_SIZE; ++j) {
                HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN[j] ==
                                   hc::modular_multiplication_prereduced_inputs(
                                                   Rmod_n[j], Rmod_n[j], n[j]));
            }
        }
    } else {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            rSquaredModN[j] = hc::modular_multiplication_prereduced_inputs(
                                                    Rmod_n[j], Rmod_n[j], n[j]);
        }
    }

    if (HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE) {
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN[j] < n[j]);
    }
    return rSquaredModN;
  }
};



template<class PTAG>
struct impl_array_get_Rsquared_mod_n<true, PTAG> {

  template <typename T, size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE static
  std::array<T, ARRAY_SIZE> call(const std::array<T, ARRAY_SIZE>& n,
                                 const std::array<T, ARRAY_SIZE>& inverse_n_modR,
                                 const std::array<T, ARRAY_SIZE>& Rmod_n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    using std::size_t;
    namespace hc = ::hurchalla;
    constexpr T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 2));
    if (HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE) {
        for (size_t j=0; j<ARRAY_SIZE; ++j) {
            HPBC_CLOCKWORK_PRECONDITION2(n[j] % 2 == 1);
            HPBC_CLOCKWORK_PRECONDITION2(n[j] > 1);
            HPBC_CLOCKWORK_PRECONDITION2(n[j] < Rdiv4);
            HPBC_CLOCKWORK_PRECONDITION2(Rmod_n[j] < n[j]);
        }
    }
    std::array<T, ARRAY_SIZE> rSquaredModN;

#ifdef HURCHALLA_TESTING_RSQUARED_MOD_N
    if (true) {
#else
    if (hc::modular_multiplication_has_slow_perf<T>()) {
#endif
        int i=0;
        std::array<T, ARRAY_SIZE> tmp(Rmod_n);   // Rmod_n ≡ 1*R (mod n)
        for (; i<4; ++i) {
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j)
                tmp[j] = hc::detail::two_times_restricted<T>::call(tmp[j], n[j]);
        }
        // at this point,  tmp ≡ 16*R (mod n)
        constexpr int bitsT = ut_numeric_limits<T>::digits;

        for (; i<bitsT/2; i*=2) {
            // use montgomery multiplication to square tmp on each iteration
            T u_hi, u_lo;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, tmp[j], tmp[j]);
                // use the same logic as MontyQuarterRange's montyREDC():
                tmp[j] = hc::REDC_incomplete(u_hi, u_lo, n[j], inverse_n_modR[j], PTAG());
                tmp[j] = static_cast<T>(tmp[j] + n[j]);
                HPBC_CLOCKWORK_ASSERT2(0 < tmp[j] && tmp[j] < static_cast<T>(2*n[j]));
            }
        }
        HPBC_CLOCKWORK_ASSERT2(i == bitsT/2);
        {
            // This final iteration was unrolled from the loop above so we can
            // use standard REDC, which will end with tmp in the range [0, n).
            T u_hi, u_lo;
            HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
                u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, tmp[j], tmp[j]);
                tmp[j] = hc::REDC_standard(u_hi, u_lo, n[j], inverse_n_modR[j], PTAG());
            }
        }
        rSquaredModN = tmp;

        if (HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE) {
            for (size_t j=0; j<ARRAY_SIZE; ++j) {
                HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN[j] ==
                                   hc::modular_multiplication_prereduced_inputs(
                                                   Rmod_n[j], Rmod_n[j], n[j]));
            }
        }
    } else {
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t j=0; j<ARRAY_SIZE; ++j) {
            rSquaredModN[j] = hc::modular_multiplication_prereduced_inputs(
                                                    Rmod_n[j], Rmod_n[j], n[j]);
        }
    }

    if (HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE) {
        for (size_t j=0; j<ARRAY_SIZE; ++j)
            HPBC_CLOCKWORK_POSTCONDITION2(rSquaredModN[j] < n[j]);
    }
    return rSquaredModN;
  }
};


}} // end namespace

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


#endif
