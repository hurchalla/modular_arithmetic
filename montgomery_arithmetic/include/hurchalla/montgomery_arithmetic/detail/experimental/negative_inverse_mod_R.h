// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_NEGATIVE_INVERSE_MODR_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_NEGATIVE_INVERSE_MODR_H_INCLUDED


#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// For discussion purposes, let R = 1<<(ut_numeric_limits<T>::digits).  For
// example if T is uint64_t, then R = 1<<64.

// minor note: uses static member functions to disallow ADL.
struct nimr_helper {
    template <int n>
    constexpr static int log2()
    {
      // PRECONDITION: n!=0 (this isn't possible to express via static_assert)
      static_assert(n>=0, "");
      static_assert(n==1 || (n/2)*2 == n, "");
      return (n<=1) ? 0 : 1 + log2<n/2>();
    }
#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
#endif

    // This algorithm is an adaptation of the algorithm described in
    // https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/low_level_api/detail/integer_inverse.pdf
    // It is an adaptation to produce the negative inverse rather than the
    // normal (positive) inverse.  The algorithm in the linked paper has to be
    // reworked from scratch (using all the same principles and the same
    // approach) to produce the negative inverse algorithm used in the function
    // below.  It's fairly straightforward to rework it for the negative inverse
    // and prove it's correct, but it's left as an exercise for the reader.
    // Note that the formula for the negative inverse of 'a' that is good to 5
    // (or 4) bits is inv_5goodbits = (3*a)^12.  Once again, you can prove the
    // correctness of this by using the same approach as the paper uses to prove
    // correctness of its formula for the *positive* inverse (good to 4 or 5
    // bits).
    //
    // Note: This alg only makes sense to use for the native integral types -
    // Newton's method becomes more efficient when larger types are required.
    template <typename T, int bits>
    static HURCHALLA_FORCE_INLINE
    typename std::enable_if<bits<=HURCHALLA_TARGET_BIT_WIDTH, T>::type
    impl_neg_inverse(T a)
    {
        static_assert(bits == ut_numeric_limits<T>::digits, "");
        static_assert(std::is_unsigned<T>::value, ""); //T native unsign integer
        HPBC_PRECONDITION2(a % 2 == 1);

        // avoid undefined behavior that could result if T is an unsigned type
        // that would be promoted to (signed) 'int'.
        using P = typename safely_promote_unsigned<T>::type;
        P b = static_cast<P>(a);
  
        P x = (3*b)^12;  // good to 5 bits, but we'll treat it as good to only 4
        static const constexpr int goodbits = 4;  // must be a power of 2
        P s = b*x;
        P y = s+1;
  
        static_assert((bits/goodbits)*goodbits == bits, "");
        static const constexpr int iterations = log2<bits/goodbits>();
        HURCHALLA_REQUEST_UNROLL_LOOP
        for (int i=0; i<iterations; ++i) {
            P t = y+1;
            y = y*y;
            x = x*t;
        }
        return static_cast<T>(x);
    }

    // This is Newton's method algorithm for the negative inverse (mod R).
    // To get the starting bits of 'x' we recurse until we can use the more
    // efficient algorithm above, at which point we switch to it.
    template <typename T, int bits>
    static HURCHALLA_FORCE_INLINE
    typename std::enable_if<!(bits<=HURCHALLA_TARGET_BIT_WIDTH), T>::type
    impl_neg_inverse(T a)
    {
        static_assert(ut_numeric_limits<T>::is_integer, "");
        static_assert(!(ut_numeric_limits<T>::is_signed), "");
        static_assert((bits/2)*2 == bits, "");
        using T2 = typename std::conditional<sized_uint<bits/2>::is_valid,
                                typename sized_uint<bits/2>::type, T>::type;
        // set x so that the lower ('bits'/2) half of the bits are good.
        T x = static_cast<T>(impl_neg_inverse<T2, bits/2>(static_cast<T2>(a)));
  
        using P = typename safely_promote_unsigned<T>::type;
        // use one step of the standard newtons method algorithm for the inverse
        // to double the number of good bits.
        return static_cast<T>(x * (2 + static_cast<P>(a)*x));
    }
};


#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4309)
#endif

// Returns the integer x satisfying  x*a ≡ -1 (mod R)
template <typename T>
T negative_inverse_mod_R(T a)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(a % 2 == 1);

    T inv= nimr_helper::impl_neg_inverse<T, ut_numeric_limits<T>::digits>(a);

    // guarantee inv*a ≡ -1 (mod R)
    using P = typename safely_promote_unsigned<T>::type;
    HPBC_POSTCONDITION2(static_cast<T>(static_cast<P>(inv) * static_cast<P>(a))
                      == static_cast<T>(static_cast<P>(0) - static_cast<P>(1)));
    return inv;
}

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


}} // end namespace

#endif
