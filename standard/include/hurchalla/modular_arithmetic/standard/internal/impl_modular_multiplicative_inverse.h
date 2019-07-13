
#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATIVE_INV_H__INCLUDED


#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


template <typename U>
U impl_modular_multiplicative_inverse(U val, U modulus)
{
    static_assert(std::is_unsigned<U>::value, "");  //U unsigned integral type
    // Precondition: modulus > 1.
    // Postcondition: Returns 0 if the inverse doesn't exist. Otherwise returns
    //    the inverse (which is never 0, given that modulus>1).


// I decided not to support modulus<=1, since it's not likely to be used and
// it complicates the return type and adds conditional branches. 
// Nevertheless I left the code here in comment, for reference:
/*
    if (modulus == 0)  // ordinarily operations modulo 0 are undefined
        return std::make_tuple((U)0, false);
    if (modulus == 1) {
        // without this "if" clause, when modulus == 1 and val == 1,
        // this function would calculate the result to be 1.  That result
        // wouldn't be completely wrong, but it isn't reduced.  We always
        // want a fully reduced result.  When modulus == 1, the fully
        // reduced result will always be 0.
        return std::make_tuple((U)0, true);
    }
*/

    using S = std::make_signed<U>::type;
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
        if (y < 0)
            return static_cast<U>(y) + modulus;
        else
            return static_cast<U>(y);
    }
    else
        return 0;
}


}}	// end namespace

#endif
