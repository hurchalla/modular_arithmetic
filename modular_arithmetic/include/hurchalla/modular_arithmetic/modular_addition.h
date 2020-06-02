
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace modular_arithmetic {


template <typename T>
T modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(modulus>0);
    HPBC_PRECONDITION(a>=0 && a<modulus);   // i.e. the input must be prereduced
    HPBC_PRECONDITION(b>=0 && b<modulus);   // i.e. the input must be prereduced

    // POSTCONDITION:
    //   Returns (a+b)%modulus.  Guarantees no overflow internally on a+b.
    
    /* We want essentially-  result = (a+b < modulus) ? a+b : a+b-modulus
       But due to potential overflow on a+b we need to write it as follows */
    T tmp = static_cast<T>(modulus - b);
    T result = (a < tmp) ? static_cast<T>(a+b) : static_cast<T>(a-tmp);
    return result;
}


}}  // end namespace

#endif
