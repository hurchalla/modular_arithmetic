
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_UNSIGNED_MULT_TO_HILO_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_UNSIGNED_MULT_TO_HILO_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/impl_unsigned_multiply_to_hilo_product.h"
#include "hurchalla/montgomery_arithmetic/internal/compiler_macros.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// Interface/contract.
// Description:  unsigned_multiply_to_hilo_product() calculates a 'double-width'
// multiplication product.  This is in contrast to a 'standard' multiply which
// drops/ignores the highest bits of the product whenever overflow occurs.
// The high-bit portion of the product is returned, and the low-bit portion is
// stored in *pLowProduct.
template <typename T>
HURCHALLA_FORCE_INLINE
T unsigned_multiply_to_hilo_product(T* pLowProduct, T a, T b)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    // POSTCONDITION: Stores the low-bits portion of the product (a*b) in
    //                *pLowProduct.
    // POSTCONDITION: Returns the high-bits portion of the product (a*b).

    return impl_unsigned_multiply_to_hilo_product(pLowProduct, a, b);
}


}} // end namespace

#endif
