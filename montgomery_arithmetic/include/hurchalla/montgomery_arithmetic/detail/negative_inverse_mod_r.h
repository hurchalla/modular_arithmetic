// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_NEGATIVE_INVERSE_MODR_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_NEGATIVE_INVERSE_MODR_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/sized_uint.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic {


namespace detail_nimr {
    template <int n>
    constexpr int log2()
    {
      // PRECONDITION: n!=0 (this isn't possible to express via static_assert)
      static_assert(n>=0, "");
      static_assert(n==1 || (n/2)*2 == n, "");
      return (n<=1) ? 0 : 1 + log2<n/2>();
    }
    #ifndef HURCHALLA_TARGET_BIT_WIDTH
    #error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
    #endif

    // This is the generalized Dumas algorithm for the negative inverse (mod R).
    // I haven't yet published my generalized form of the Dumas algorithm, but
    // the Dumas algorithm comes from  https://arxiv.org/abs/1209.6626
    // The closest information available at the moment is from Marc Reynolds at
    // http://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html
    // However, Reynolds presents a straightforward adaptation of Dumas's
    // algorithm.  This generalized form is a slightly different algo.
    //
    // Note: Dumas's alg only makes sense to use for the native integral types -
    // Newton's method becomes more efficient when larger types are required.
    template <typename T, int bits>
    HURCHALLA_FORCE_INLINE
    typename std::enable_if<bits<=HURCHALLA_TARGET_BIT_WIDTH, T>::type
    impl_neg_inverse(T a)
    {
        namespace ma = hurchalla::modular_arithmetic;
        static_assert(bits == ma::ma_numeric_limits<T>::digits, "");
        static_assert(std::is_unsigned<T>::value, ""); //T native unsign integer
        HPBC_PRECONDITION2(a % 2 == 1);
        HPBC_PRECONDITION2(a > 1);

        // avoid undefined behavior that could result if T is an unsigned type
        // that would be promoted to (signed) 'int'.
        using U = typename safely_promote_unsigned<T>::type;
        U b = static_cast<U>(a);
  
        U x = (3*b)^12;  // good to 5 bits, but we'll treat it as good to only 4
        static const constexpr int goodbits = 4;  // must be a power of 2
        U s = b*x;
        U y = s+1;
  
        static_assert((bits/goodbits)*goodbits == bits, "");
        static const constexpr int iterations = log2<bits/goodbits>();
        HURCHALLA_REQUEST_UNROLL_LOOP
        for (int i=0; i<iterations; ++i) {
            U t = y+1;
            y = y*y;
            x = x*t;
        }
        return static_cast<T>(x);
    }

    // This is Newton's method algorithm for the negative inverse (mod R).
    // To get the starting bits of 'x' we recurse until we use Dumas's method
    // (it's more efficient than Newton's method for native integer types).
    template <typename T, int bits>
    HURCHALLA_FORCE_INLINE
    typename std::enable_if<!(bits<=HURCHALLA_TARGET_BIT_WIDTH), T>::type
    impl_neg_inverse(T a)
    {
        static_assert((bits/2)*2 == bits, "");
        using T2 = typename sized_uint<bits/2>::type;
        using T3 = typename std::conditional<!(std::is_same<T2,void>::value),
                                         T2, T>::type;
        // set x so that the lower ('bits'/2) half of the bits are good.
        T x = static_cast<T>(impl_neg_inverse<T3, bits/2>(static_cast<T3>(a)));
  
        // use one step of the standard newtons method algorithm for the inverse
        // to double the number of good bits.
        return x*(static_cast<T>(2) + a*x);
    }
}


// Returns the integer x satisfying  x*a ≡ -1 (mod R)
template <typename T>
T negative_inverse_mod_r(T a)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(a % 2 == 1);
    HPBC_PRECONDITION2(a > 1);

    T inv= detail_nimr::impl_neg_inverse<T,ma::ma_numeric_limits<T>::digits>(a);

    // guarantee inv*a ≡ -1 (mod R)
    using U = typename safely_promote_unsigned<T>::type;
    HPBC_POSTCONDITION2(static_cast<T>(static_cast<U>(inv) * static_cast<U>(a))
                      == static_cast<T>(static_cast<U>(0) - static_cast<U>(1)));
    return inv;
}


}} // end namespace

#endif
