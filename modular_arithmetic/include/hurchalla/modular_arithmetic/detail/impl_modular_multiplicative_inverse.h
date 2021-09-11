// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H_INCLUDED


#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


template <typename T> HURCHALLA_FORCE_INLINE
T impl_modular_multiplicative_inverse(T val, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    // I decided not to support modulus<=1, since it's not likely to be used and
    // it complicates the return type and adds conditional branches.
    HPBC_PRECONDITION2(modulus > 1);

    // POSTCONDITION: Returns 0 if the inverse doesn't exist. Otherwise returns
    //    the inverse (which is never 0, given that modulus>1).

    using U = typename safely_promote_unsigned<T>::type;
    using S = typename extensible_make_signed<U>::type;

    // The following algorithm is adapted from Figure 6 of
    // https://jeffhurchalla.com/2018/10/13/implementing-the-extended-euclidean-algorithm-with-unsigned-inputs/
    // calculating only what is needed for the modular multiplicative inverse.
    S y1=0;
    U a1=modulus;
    {
        S y0=1;
        U a2=val;
        U q=0;
        while (a2 != 0) {
            S y2 = static_cast<S>(y0 - static_cast<S>(q)*y1);
            y0=y1;
            y1=y2;
            U a0=a1;
            a1=a2;

            q = static_cast<U>(a0/a1);
            a2 = static_cast<U>(a0 - q*a1);
        }
    }
    S y = y1;
    U gcd = a1;
    if (gcd == 1) {
#if 0
        U inv = (y < 0) ? static_cast<U>(static_cast<U>(y) + modulus) :
                                                              static_cast<U>(y);
#else
        U inv = static_cast<U>(y);
        HURCHALLA_CMOV(y < 0, inv, static_cast<U>(static_cast<U>(y) + modulus));
#endif
        HPBC_POSTCONDITION2(inv < modulus);
        return static_cast<T>(inv);
    } else
        return 0;
}


}}  // end namespace

#endif
