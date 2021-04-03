// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_INVERSE_MOD_R_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_INVERSE_MOD_R_H_INCLUDED


#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// For discussion purposes, let the unlimited precision constant R equal
// 2^(ut_numeric_limits<T>::digits).  For example when T is uint64_t, R = 2^64.


namespace detail_invmR {
    template <int n>
    constexpr int log2()
    {
        // PRECONDITION: n!=0 (this isn't possible to express via static_assert)
        static_assert(n>=0, "");
        static_assert(n==1 || (n/2)*2 == n, "");
        return (n<=1) ? 0 : 1 + log2<n/2>();
    }
}  // end namespace detail_invmR

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
#endif

// This is the generalized Dumas algorithm for the inverse (mod R).
// TODO:
// I haven't yet published my generalized form of the Dumas algorithm, but the
// Dumas algorithm comes from  https://arxiv.org/abs/1209.6626
// The closest information available at the moment is from Marc Reynolds at
// http://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html
// However, Reynolds presents a straightforward adaptation of Dumas's algorithm.
// This generalized form is a slightly different algo.
//
// Note: Dumas's alg only makes sense to use for the native integral types -
// Newton's method becomes more efficient when larger types are required.
template <typename T, int bits>
HURCHALLA_FORCE_INLINE
typename std::enable_if<(bits <= HURCHALLA_TARGET_BIT_WIDTH), T>::type
impl_inverse_mod_R(T a)
{
    static_assert(bits == ut_numeric_limits<T>::digits, "");
    static_assert(std::is_unsigned<T>::value, "");  // T native unsigned integer
    HPBC_PRECONDITION2(a % 2 == 1);

    // avoid undefined behavior that could result if T is an unsigned type that
    // would be promoted to (signed) 'int'.
    using P = typename safely_promote_unsigned<T>::type;
    P b = static_cast<P>(a);

    P x = (3*b)^2;  // good to 5 bits, but we'll treat it as good to only 4
    static const constexpr int goodbits = 4;  // must be a power of 2
    P s = b*x;
    P y = 1-s;

    static_assert((bits/goodbits)*goodbits == bits, "");
    static const constexpr int iterations = detail_invmR::log2<bits/goodbits>();
    HURCHALLA_REQUEST_UNROLL_LOOP
    for (int i=0; i<iterations; ++i) {
        P t = y+1;
        y = y*y;
        x = x*t;
    }
    return static_cast<T>(x);
}

// This is Newton's method algorithm for the inverse (mod R).
// To get the starting bits of 'x' we recurse until we use Dumas's method (it's
// more efficient than Newton's method for native integer types).
template <typename T, int bits>
HURCHALLA_FORCE_INLINE
typename std::enable_if<!(bits <= HURCHALLA_TARGET_BIT_WIDTH), T>::type
impl_inverse_mod_R(T a)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert((bits/2)*2 == bits, "");
    using T2 = typename std::conditional<sized_uint<bits/2>::is_valid,
                                    typename sized_uint<bits/2>::type, T>::type;
    HPBC_PRECONDITION2(a % 2 == 1);

    // set x so that the lower ('bits'/2) half of the bits are good.
    T x = static_cast<T>(impl_inverse_mod_R<T2, bits/2>(static_cast<T2>(a)));

    using P = typename safely_promote_unsigned<T>::type;
    // use one step of the standard newtons method algorithm for the inverse to
    // double the number of good bits.
    return static_cast<T>(x * (2 - static_cast<P>(a)*x));
}


}} // end namespace

#endif
