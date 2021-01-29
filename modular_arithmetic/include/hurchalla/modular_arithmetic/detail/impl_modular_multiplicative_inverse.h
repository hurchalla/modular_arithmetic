// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H_INCLUDED


#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// Enabled only for signed (integral) types T.
template <typename T> HURCHALLA_FORCE_INLINE
typename std::enable_if<ut_numeric_limits<T>::is_signed, T>::type
impl_modular_multiplicative_inverse(T val, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(val >= 0);
    HPBC_PRECONDITION2(modulus > 1);

    // POSTCONDITION: Returns 0 if the inverse doesn't exist. Otherwise returns
    //    the inverse (which is never 0, given that modulus>1).


// I decided not to support modulus<=1, since it's not likely to be used and
// it complicates the return type and adds conditional branches. 
// Nevertheless I left the code here in comment, for reference:
/*
    if (modulus == 0)  // ordinarily operations modulo 0 are undefined
        return std::make_tuple(static_cast<T>(0), false);
    if (modulus == 1) {
        // without this "if" clause, when modulus == 1 and val == 1,
        // this function would calculate the result to be 1.  That result
        // wouldn't be completely wrong, but it isn't reduced.  We always
        // want a fully reduced result.  When modulus == 1, the fully
        // reduced result will always be 0.
        return std::make_tuple(static_cast<T>(0), true);
    }
*/

    // The following algorithm is simplified from Figure 1 of
    // https://jeffhurchalla.com/2018/10/13/implementing-the-extended-euclidean-algorithm-with-unsigned-inputs/
    // calculating only what is needed for the modular multiplicative inverse.
    T gcd, y;
    {
       T y0=0, a0=modulus;
       T y1=1, a1=val;
    
       while (a1 != 0) {
          T q = static_cast<T>(a0/a1);
          T a2 = static_cast<T>(a0 - q*a1);
          T y2 = static_cast<T>(y0 - q*y1);
          y0=y1; a0=a1;
          y1=y2; a1=a2;
       }
       y = y0;
       gcd = a0;
    }

    if (gcd == 1) {
        T inv = (y < 0) ? static_cast<T>(y + modulus) : y;
        HPBC_POSTCONDITION2(inv>=0 && inv<modulus);
        return inv;
    } else
        return 0;
}


// Enabled only for unsigned (integral) types U.
template <typename U> HURCHALLA_FORCE_INLINE
typename std::enable_if<!(ut_numeric_limits<U>::is_signed), U>::type
impl_modular_multiplicative_inverse(U val, U modulus)
{
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");
    HPBC_PRECONDITION2(modulus>1);

    // POSTCONDITION: Returns 0 if the inverse doesn't exist. Otherwise returns
    //    the inverse (which is never 0, given that modulus>1).


// I decided not to support modulus<=1, since it's not likely to be used and
// it complicates the return type and adds conditional branches. 
// Nevertheless I left the code here in comment, for reference:
/*
    if (modulus == 0)  // ordinarily operations modulo 0 are undefined
        return std::make_tuple(static_cast<T>(0), false);
    if (modulus == 1) {
        // without this "if" clause, when modulus == 1 and val == 1,
        // this function would calculate the result to be 1.  That result
        // wouldn't be completely wrong, but it isn't reduced.  We always
        // want a fully reduced result.  When modulus == 1, the fully
        // reduced result will always be 0.
        return std::make_tuple(static_cast<T>(0), true);
    }
*/
    using S = typename extensible_make_signed<U>::type;

    // The following algorithm is adapted from Figure 6 of
    // https://jeffhurchalla.com/2018/10/13/implementing-the-extended-euclidean-algorithm-with-unsigned-inputs/
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
        U inv = (y < 0) ? static_cast<U>(static_cast<U>(y) + modulus) :
                                                              static_cast<U>(y);
        HPBC_POSTCONDITION2(inv<modulus);
        return inv;
    } else
        return 0;
}


}}  // end namespace

#endif
