
#ifndef HURCHALLA_MODULAR_ARITHMETIC_MONTGOMERY_FORM_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_MONTGOMERY_FORM_H__INCLUDED


#include "hurchalla/modular_arithmetic/montgomery/internal/montgomerydefault.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace modular_montgomery {


// When using the default MontyType, T must be signed or unsigned integral type.
// A custom MontyType (SIMD perhaps) might require a different type T.
template<typename T, class MontyType = MontgomeryDefault<T>::type>
class MontgomeryForm final {
    MontyType impl;
    using U = MontyType::value_type;
    static_assert(sizeof(U) >= sizeof(T));
public:
    using V = MontyType::montvalue_type; // MontgomeryValue<U>;

    explicit MontgomeryForm(T modulus) noexcept : impl(static_cast<U>(modulus))
    {
        precondition(modulus & 1 == 1);  // modulus must be odd
        precondition(modulus > 1);
    }
    MontgomeryForm(const MontgomeryForm&) = delete;
    MontgomeryForm& operator=(const MontgomeryForm&) = delete;

    // Returns the converted value of the standard number 'a' into monty form.
    // Requires 0 <= a < modulus.  The return value might not be canonical -
    // call getCanonicalForm() if you need to use it in comparisons.
    V convertIn(T a) const noexcept
    {
        precondition(0 <= a);
        precondition(a < static_cast<T>(impl.getModulus()));
        return impl.convertIn(static_cast<U>(a));
    }
    // Converts (montgomery value) x into a "normal" number; returns the result.
    // Guarantees 0 <= result < modulus.
    T convertOut(V x) const noexcept
    {
        T a = static_cast<T>(impl.convertOut(x));
        postcondition(0 <= a);
        postcondition(a < static_cast<T>(impl.getModulus()));
        return a;
    }

    // Returns a unique (canonical) value representing the equivalence class of
    // x modulo the modulus.  Note that this return value can be used to test
    // for equality with another canonical value in montgomery form.
    V getCanonicalForm(V x) const noexcept { return impl.getCanonicalForm(x); }

    // Returns the canonical converted value of 1 in montgomery form.
    V getUnityValue() const noexcept
    {
        V ret = impl.getUnityValue();
        postcondition(impl.isCanonical(ret));
        return ret;
    }
    // Returns the canonical converted value of 0 in montgomery form.
    V getZeroValue() const noexcept
    {
        V ret = impl.getZeroValue();
        postcondition(impl.isCanonical(ret));
        return ret;
    }
    // Returns the canonical converted value of modulus-1 (or -1) in monty form.
    V getNegativeOneValue() const noexcept
    {
        V ret = impl.getNegativeOneValue();
        postcondition(impl.isCanonical(ret));
        return ret;
    }

    // Returns the modular product of (the montgomery values) x and y.  The
    // return value is in montgomery form but might not be canonical - call
    // getCanonicalForm() to use it in comparisons.
//    FORCE_INLINE V inline_multiply(V x, V y) const noexcept { 
//        return impl.inline_multiply(x, y); 
//    }
    // Same discussion/return as inline_multiply(). This function is not inline.
    V multiply(V x, V y) const noexcept { return impl.multiply(x, y); }

    // Returns the modular sum of (the montgomery values) x and y.  The return
    // value is in montgomery form but might not be canonical - call
    // getCanonicalForm() to use it in comparisons.
    V add(V x, V y) const noexcept { return impl.add(x, y); }

    // Returns the modular difference of (the montgomery values) x and y.  More
    // precisely, x minus y.  The return value is in montgomery form but might
    // not be canonical - call getCanonicalForm() to use it in comparisons.
    V subtract(V x, V y) const noexcept { return impl.subtract(x, y); }
};


}} // end namespace

#endif
