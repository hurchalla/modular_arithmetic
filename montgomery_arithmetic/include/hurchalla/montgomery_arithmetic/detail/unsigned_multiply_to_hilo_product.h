// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_UNSIGNED_MULT_TO_HILO_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_UNSIGNED_MULT_TO_HILO_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/impl_unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"

namespace hurchalla { namespace montgomery_arithmetic {


// unsigned_multiply_to_hilo_product() calculates a 'double-width'
// multiplication product.  This behavior differs from a 'standard' multiply
// which drops/ignores the highest bits of the product whenever overflow occurs.
//
// Returns the high-bit portion of the product, and stores the low-bit portion
// in *pLowProduct.
template <typename T>
HURCHALLA_FORCE_INLINE
T unsigned_multiply_to_hilo_product(T* pLowProduct, T a, T b)
{
    namespace ut = hurchalla::util;
    static_assert(ut::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut::ut_numeric_limits<T>::is_signed), "");
    // POSTCONDITION: Stores the low-bits portion of the product (a*b) in
    //                *pLowProduct.
    // POSTCONDITION: Returns the high-bits portion of the product (a*b).

    return impl_unsigned_multiply_to_hilo_product(pLowProduct, a, b);
}


}} // end namespace

#endif
