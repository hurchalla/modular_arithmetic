
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_ADDITION_H__INCLUDED


#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace modular_arithmetic {


template <typename T>
T modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    precondition(modulus>0);
    precondition(a>=0 && a<modulus);   // i.e. the input must be prereduced
    precondition(b>=0 && b<modulus);   // i.e. the input must be prereduced
    // Postcondition:
    //   Returns (a+b)%modulus.  Guarantees no overflow internally on a+b.
    
    /* We want essentially-  result = (a+b < modulus) ? a+b : a+b-modulus
       But due to potential overflow on a+b we need to write it as follows */
    T tmp = modulus - b;
    T result = (a < tmp) ? a+b : a-tmp;
    return result;
}


}}  // end namespace

#endif
