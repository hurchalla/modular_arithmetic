
#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H_INCLUDED


#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>
#include <limits>

namespace hurchalla { namespace modular_arithmetic {


// Enabled only for signed (integral) types T.
template <typename T>
typename std::enable_if<std::numeric_limits<T>::is_signed, T>::type
impl_modular_multiplicative_inverse(T val, T modulus)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(std::numeric_limits<T>::is_signed, "");
    precondition(val >= 0);
    precondition(modulus > 1);
    // Postcondition: Returns 0 if the inverse doesn't exist. Otherwise returns
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
          T q = a0/a1;
          T a2 = a0 - q*a1;
          T y2 = y0 - q*y1;
          y0=y1; a0=a1;
          y1=y2; a1=a2;
       }
       y = y0;
       gcd = a0;
    }

    if (gcd == 1) {
        T inv = (y < 0) ? y + modulus : y;
        postcondition(inv>=0 && inv<modulus);
        return inv;
    }
    else
        return 0;
}


// Enabled only for unsigned (integral) types U.
template <typename U>
typename std::enable_if<!(std::numeric_limits<U>::is_signed), U>::type
impl_modular_multiplicative_inverse(U val, U modulus)
{
    static_assert(std::numeric_limits<U>::is_integer, "");
    static_assert(!(std::numeric_limits<U>::is_signed), "");
    precondition(modulus>1);
    // Postcondition: Returns 0 if the inverse doesn't exist. Otherwise returns
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

    // std::make_signed only compiles for types known to std::type_traits, which
    // means basically the native integral types and nothing else.
    // If you have a compile error here, you can probably create a non-template
    // overload of this function that is specific to the unsigned integral type
    // that you wanted to use.  The body of the overload should essentially be a
    // copy of this function, but with the line below changed to explicitly
    // specify the signed version of your type (instead of std::make_signed).
    using S = typename std::make_signed<U>::type;

    // The following algorithm is adapted from Figure 6 of
    // https://jeffhurchalla.com/2018/10/13/implementing-the-extended-euclidean-algorithm-with-unsigned-inputs/
    S y1=0;
    U a1=modulus;
    {
        S y0=1;
        U a2=val;
        U q=0;
        while (a2 != 0) {
            S y2 = y0 - static_cast<S>(q)*y1;
            y0=y1;
            y1=y2;
            U a0=a1;
            a1=a2;

            q = a0/a1;
            a2 = a0 - q*a1;
        }
    }
    S y = y1;
    U gcd = a1;
    if (gcd == 1) {
        U inv = (y < 0) ? static_cast<U>(y) + modulus : static_cast<U>(y);
        postcondition(inv<modulus);
        return inv;
    }
    else
        return 0;
}


}}  // end namespace

#endif
