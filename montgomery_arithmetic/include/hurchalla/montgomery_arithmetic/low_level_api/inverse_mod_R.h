// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_INVERSE_MOD_R_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_INVERSE_MOD_R_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/detail/impl_inverse_mod_R.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// For discussion purposes, let the unlimited precision constant R equal
// 2^(ut_numeric_limits<T>::digits).  For example when T is uint64_t, R = 2^64.

// Returns the integer x satisfying  x*a ≡ 1 (mod R)
template <typename T>
T inverse_mod_R(T a)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(a % 2 == 1);

    T inv = detail::impl_inverse_mod_R<T, ut_numeric_limits<T>::digits>(a);

    // guarantee inv*a ≡ 1 (mod R)
    using P = typename safely_promote_unsigned<T>::type;
    HPBC_POSTCONDITION2(static_cast<T>(static_cast<P>(inv) * static_cast<P>(a))
                        == static_cast<T>(1));
    return inv;
}


} // end namespace

#endif
