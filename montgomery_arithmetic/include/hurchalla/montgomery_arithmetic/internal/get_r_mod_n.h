
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_GET_R_MOD_N_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_GET_R_MOD_N_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// For discussion purposes, given an unsigned integral type T, let
// R = 2^(std::numeric_limits<T>::digits). For example: if T is uint64_t
// then R = 2^64.
//
// Returns r_mod_n == R%n
template <typename T>
T get_r_mod_n(T n)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(n > 1);

    // Assign a tmp T variable rather than directly using the intermediate
    // expression, in order to avoid a negative value (and a wrong answer) in
    // cases where 'n' would be promoted to type 'int'
    T tmp = static_cast<T>(0) - n;

    // Compute R%n.  For example, if R==2^64, arithmetic wraparound behavior of
    // the unsigned integral type T results in (0 - n) representing (2^64 - n).
    // Thus, r_mod_n = R%n == (2^64)%n == (2^64 - n)%n == (0-n)%n
    T r_mod_n = tmp % n;

    // Since n is odd and > 1, n does not divide R==2^x.  Thus, r_mod_n != 0
    HPBC_POSTCONDITION2(0 < r_mod_n && r_mod_n < n);
    return r_mod_n;
}


}} // end namespace

#endif
