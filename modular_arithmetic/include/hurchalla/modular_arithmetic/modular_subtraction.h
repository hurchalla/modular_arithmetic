
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MODULAR_SUBTRACTION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MODULAR_SUBTRACTION_H_INCLUDED


#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace modular_arithmetic {


template <typename T>
T modular_subtraction_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    precondition(modulus>0);
    precondition(a>=0 && a<modulus);   // i.e. the input must be prereduced
    precondition(b>=0 && b<modulus);   // i.e. the input must be prereduced
    // Postcondition:
    //   Returns (a-b)%modulus.  Guarantees no underflow internally on a-b.

    /* We want essentially-  result = (a-b < 0) ? a-b+modulus : a-b
        But due to potential overflow on a-b we need to write it as follows */
    T result = (a>=b) ? a-b : modulus-(b-a);
    return result;

    /* The branches for both modular_subtraction and modular_addition are
        usually unpredictable, so the following *might* be faster (it
        avoids a conditional branch)
        T maybe_addend = modulus & (static_cast<T>(0) - static_cast<T>(a < b));
        T result = (a - b) + maybe_addend;
        This is probably worse than a conditional move (cmov) though.
        Really, what we would like is a way to designate that a branch is
        unpredictable, like clang's __builtin_unpredictable().  Presumably the
        compiler would use conditional moves in that case.
    */
}


}}  // end namespace

#endif
