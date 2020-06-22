// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_SAFELY_PROMOTE_UNSIGNED_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_SAFELY_PROMOTE_UNSIGNED_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic {


// safely_promote_unsigned<T> is intended to protect against the undefined
// behavior and unexpected results that can arise from unsigned integral
// promotion in C++.  For details on these issues, see
// https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/

// If an unsigned type T would get promoted to (signed) 'int', we want to make
// sure that the type safely_promote_unsigned<T> provides is 'unsigned int'.
// Otherwise the type it provides is T.


// Note that this primary template will be enabled for non-native unsigned
// integral types.  C++ never promotes non-native integer types, so this primary
// template just provides back the type T.
template <typename T, typename Enable = void>
struct safely_promote_unsigned {
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    using type = T;
};

// Note that this specialization will be enabled for native unsigned integral
// types (std::is_unsigned<T>::value is true only for native unsigned int types)
template <typename T>
struct safely_promote_unsigned<T, 
                    typename std::enable_if<std::is_unsigned<T>::value>::type> {
    using type = typename std::make_unsigned<
                           decltype(static_cast<T>(1)*static_cast<T>(1))>::type;
};


}} // end namespace

#endif
