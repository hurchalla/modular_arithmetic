
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MAKE_SAFE_UNSIGNED_INTEGER_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MAKE_SAFE_UNSIGNED_INTEGER_H_INCLUDED


#include <type_traits>
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// make_safe_unsigned_integer<T> is intended to protect against the undefined
// behavior and unexpected results that can arise from unsigned integral
// promotion in C++.  For details on these issues, see
// https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/


// If an unsigned type T would get promoted to (signed) 'int', we want to make
// sure that the type make_safe_unsigned_integer<T> provides is 'unsigned int'.
// Otherwise the type it provides is T.


// Note that this primary template will be enabled for non-native unsigned
// integral types.  C++ never promotes non-native integer types, so this primary
// template just provides back the type T.
template <typename T, typename Enable = void>
struct make_safe_unsigned_integer {
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    using type = T;
};

// Note that this specialization will be enabled for native unsigned integral
// types (std::is_unsigned<T>::value is true only for native unsigned int types)
template <typename T>
struct make_safe_unsigned_integer<T, 
                    typename std::enable_if<std::is_unsigned<T>::value>::type> {
    using type = typename std::make_unsigned<decltype((T)1*(T)1)>::type;
};


}} // end namespace

#endif
