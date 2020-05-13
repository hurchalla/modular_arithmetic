
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/MontgomeryDefault.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// When using the default MontyType, T must be signed or unsigned integral type.
// A custom MontyType(SIMD perhaps) might have different requirements for type T
template<typename T, class MontyType = typename MontgomeryDefault<T>::type>
class MontgomeryForm final {
    MontyType impl;
    using U = typename MontyType::template_param_type;
    static_assert(std::numeric_limits<U>::is_integer, "");
    static_assert(!(std::numeric_limits<U>::is_signed), "");
    static_assert(std::numeric_limits<U>::digits >=
                  std::numeric_limits<T>::digits, "");
public:
    using V = typename MontyType::montvalue_type; // MontgomeryValue<U>;

    explicit MontgomeryForm(T modulus) : impl(static_cast<U>(modulus))
    {
        precondition(modulus % 2 == 1);  // modulus must be odd
        precondition(modulus > 1);
    }
    MontgomeryForm(const MontgomeryForm&) = delete;
    MontgomeryForm& operator=(const MontgomeryForm&) = delete;

    // Returns the converted value of the standard number 'a' into monty form.
    // Requires 0 <= a < modulus.  The return value might not be canonical -
    // call getCanonicalForm() if you need to use it in comparisons.
    V convertIn(T a) const
    {
        precondition(a >= 0);
        precondition(a < static_cast<T>(impl.getModulus()));
        return impl.convertIn(static_cast<U>(a));
    }
    // Converts (montgomery value) x into a "normal" number; returns the result.
    // Guarantees 0 <= result < modulus.
    T convertOut(V x) const
    {
        T a = static_cast<T>(impl.convertOut(x));
        postcondition(a >= 0);
        postcondition(a < static_cast<T>(impl.getModulus()));
        return a;
    }

    // Returns a unique (canonical) value representing the equivalence class of
    // x modulo the modulus.  Note that this return value can be used to test
    // for equality with another canonical value in montgomery form.
    V getCanonicalForm(V x) const { return impl.getCanonicalForm(x); }

    // Returns the canonical converted value of 1 in montgomery form.
    V getUnityValue() const
    {
        V ret = impl.getUnityValue();
        postcondition(impl.isCanonical(ret));
        return ret;
    }
    // Returns the canonical converted value of 0 in montgomery form.
    V getZeroValue() const
    {
        V ret = impl.getZeroValue();
        postcondition(impl.isCanonical(ret));
        return ret;
    }
    // Returns the canonical converted value of modulus-1 (or -1) in monty form.
    V getNegativeOneValue() const
    {
        V ret = impl.getNegativeOneValue();
        postcondition(impl.isCanonical(ret));
        return ret;
    }

    // Returns the modular product of (the montgomery values) x and y.  The
    // return value is a montgomery value but might not be canonical - call
    // getCanonicalForm() to use it in comparisons.
    V multiply(V x, V y) const { return impl.multiply(x, y); }

//    FORCE_INLINE V inline_multiply(V x, V y) const { 
//        return impl.multiply(x, y); 
//    }

    // Returns the modular product of (the montgomery value) x multiplied by x.
    // The return value is a montgomery value but might not be canonical - call
    // getCanonicalForm() to use it in comparisons.
    V square(V x) const { return impl.multiply(x, x); }

    // Calculates and returns the modular exponentiation of the montgomery value
    // 'base' to the power of (the type T variable) 'exponent'.  The return
    // value is a mongomery value but might not be canonical - call
    // getCanonicalForm() to use it in comparisons.
    V pow(V base, T exponent)
    {
        // This is a slightly optimized version of Algorithm 14.76, from
        // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
        V result;
        if (exponent & static_cast<T>(1))
            result = base;
        else
            result = impl.getUnityValue();
        while (exponent > static_cast<T>(1))
        {
            exponent = exponent >> static_cast<T>(1);
            base = impl.multiply(base, base);
            if (exponent & static_cast<T>(1))
                result = impl.multiply(result, base);
        }
        return result;
    }

    // Returns the modular sum of (the montgomery values) x and y.  The return
    // value might not be canonical- use getCanonicalForm() for comparisons.
    V add(V x, V y) const { return impl.add(x, y); }

    // Returns the modular difference of (the montgomery values) x and y.  More
    // precisely, x minus y.  The return value might not be canonical- use
    // getCanonicalForm() for comparisons.
    V subtract(V x, V y) const { return impl.subtract(x, y); }

//    bool isValid(V x) const { return impl.isValid(x); }
//    bool isReduced(T a) const { return impl.isReduced(a); }
//    bool isCanonical(V x) const { return impl.isCanonical(x); }
//    T getModulus() const { return impl.getModulus(x); }
};


}} // end namespace

#endif
