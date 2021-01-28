// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_GET_RSQUARED_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_GET_RSQUARED_MOD_N_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// For discussion purposes, let the unlimited precision constant R represent
// R = 2^(ut_numeric_limits<T>::digits).  For example, if T is uint64_t, then
// R = 2^64.

// Compute (R*R) % n
template <typename T, class MTAG = FullrangeTag> HURCHALLA_FORCE_INLINE
T impl_get_Rsquared_mod_n(T n, T inverse_n_modR, T Rmod_n, MTAG = MTAG())
{
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error HURCHALLA_TARGET_BIT_WIDTH must be defined
#endif
#if defined(HURCHALLA_TARGET_ISA_X86_32) || defined(HURCHALLA_TARGET_ISA_X86_64)
    // x86 has a division instruction that has a dividend parameter that is
    // twice the CPU word size (word size == HURCHALLA_TARGET_BIT_WIDTH).
    // Modular multiplication produces a temporary product that is twice its
    // operand bit width, and also divides a temporary dividend that is twice
    // the operand bit width.  That's why we flag if we're on x86.
    constexpr bool is_x86 = true;
#else
    constexpr bool is_x86 = false;
#endif
#ifdef HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE
    constexpr bool no_native_divide = true;
#else
    constexpr bool no_native_divide = false;
#endif
#ifdef HURCHALLA_TESTING_RSQUARED_MOD_N
    constexpr bool is_a_test = true;
#else
    constexpr bool is_a_test = false;
#endif
    constexpr int bitsT = ut_numeric_limits<T>::digits;

    T rSquaredModN;
    if (no_native_divide || (bitsT > HURCHALLA_TARGET_BIT_WIDTH) ||
              ((bitsT == HURCHALLA_TARGET_BIT_WIDTH) && !is_x86) || is_a_test) {
        HPBC_ASSERT2(Rmod_n < n);
        T tmp = Rmod_n;   // Rmod_n == 1*R (mod n)
        int i=0;
        for (; i<4; ++i)
            tmp = modular_addition_prereduced_inputs(tmp, tmp, n);
        // at this point,  tmp == 16*R (mod n)
        for (; i<bitsT; i*=2) {
            // use montgomery multiplication to square tmp on each iteration
            T u_hi, u_lo;
            u_hi = unsigned_multiply_to_hilo_product(&u_lo, tmp, tmp);
            tmp = REDC(u_hi, u_lo, n, inverse_n_modR, MTAG(), LowlatencyTag());
        }
        HPBC_ASSERT2(i == bitsT);
        // we should now have  tmp ≡ R*R (mod n).

        // REDC's postcondition guarantees the following
        HPBC_ASSERT2((std::is_same<MTAG, FullrangeTag>::value ||
                      std::is_same<MTAG, HalfrangeTag>::value) ?
                     tmp < n : 0 < tmp && tmp < 2*n);
        if (!(std::is_same<MTAG, FullrangeTag>::value ||
                      std::is_same<MTAG, HalfrangeTag>::value)) {
            if (tmp >= n)
                tmp = static_cast<T>(tmp - n);  // fully reduce tmp, mod n.
        }
        rSquaredModN = tmp;
        HPBC_POSTCONDITION2(rSquaredModN ==
                   modular_multiplication_prereduced_inputs(Rmod_n, Rmod_n, n));
    } else {
        rSquaredModN = modular_multiplication_prereduced_inputs(
                                                             Rmod_n, Rmod_n, n);
    }

    HPBC_POSTCONDITION2(rSquaredModN < n);
    return rSquaredModN;
}


}} // end namespace

#endif
