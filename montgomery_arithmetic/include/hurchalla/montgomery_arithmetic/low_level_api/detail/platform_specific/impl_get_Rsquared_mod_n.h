// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_GET_RSQUARED_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_GET_RSQUARED_MOD_N_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace detail {


// For discussion purposes, let the unlimited precision constant R represent
// R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t, then
// R = 1<<64.

// Compute (R*R) % n
// Minor note: uses static member function to disallow ADL.
struct impl_get_Rsquared_mod_n {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T n, T inverse_n_modR, T Rmod_n)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    namespace hc = ::hurchalla;
    T rSquaredModN;
#ifdef HURCHALLA_TESTING_RSQUARED_MOD_N
    if (true) {
#else
    if (hc::modular_multiplication_has_slow_perf<T>()) {
#endif
        HPBC_ASSERT2(Rmod_n < n);
        T tmp = Rmod_n;   // Rmod_n ≡ 1*R (mod n)
        int i=0;
        for (; i<4; ++i)
            tmp = hc::modular_addition_prereduced_inputs(tmp, tmp, n);
        // at this point,  tmp ≡ 16*R (mod n)
        constexpr int bitsT = ut_numeric_limits<T>::digits;
        for (; i<bitsT; i*=2) {
            // use montgomery multiplication to square tmp on each iteration
            T u_hi, u_lo;
            u_hi = hc::unsigned_multiply_to_hilo_product(u_lo, tmp, tmp);
            tmp = hc::REDC_standard(
                                u_hi, u_lo, n, inverse_n_modR, LowlatencyTag());
        }
        HPBC_ASSERT2(i == bitsT);
        // We should now have  tmp ≡ R*R (mod n).
        // REDC_standard's postcondition guarantees the following:
        HPBC_ASSERT2(tmp < n);

        rSquaredModN = tmp;
        HPBC_POSTCONDITION2(rSquaredModN ==
               hc::modular_multiplication_prereduced_inputs(Rmod_n, Rmod_n, n));
    } else {
        rSquaredModN = hc::modular_multiplication_prereduced_inputs(
                                                             Rmod_n, Rmod_n, n);
    }

    HPBC_POSTCONDITION2(rSquaredModN < n);
    return rSquaredModN;
  }
};


}} // end namespace

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


#endif
