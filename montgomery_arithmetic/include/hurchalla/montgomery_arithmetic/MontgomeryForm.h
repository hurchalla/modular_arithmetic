// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


// When using the default MontyType, T must be signed or unsigned integral type.
// A custom MontyType may have different requirements for type T (e.g. that T is
// an unsigned integral type).
template<typename T, class MontyType = typename MontgomeryDefault<T>::type>
class MontgomeryForm final {
    const MontyType impl;
    using U = typename MontyType::template_param_type;
    static_assert(modular_arithmetic::ma_numeric_limits<U>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<U>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<U>::digits >=
                  modular_arithmetic::ma_numeric_limits<T>::digits, "");
public:
    using MontgomeryValue = typename MontyType::montvalue_type;
    using T_type = T;
    class CanonicalValue : public MontgomeryValue {
        friend MontgomeryForm;
        explicit CanonicalValue(MontgomeryValue val) : MontgomeryValue(val) {}
    public:
        CanonicalValue() : MontgomeryValue() {}
        friend bool operator==(const CanonicalValue& x, const CanonicalValue& y)
        {
            return x.value == y.value;
        }
        friend bool operator!=(const CanonicalValue& x, const CanonicalValue& y)
        {
            return !(x == y);
        }
    };

    explicit MontgomeryForm(T modulus) : impl(static_cast<U>(modulus))
    {
        HPBC_PRECONDITION(modulus % 2 == 1);  // modulus must be odd
        HPBC_PRECONDITION(modulus > 1);
    }
    MontgomeryForm(const MontgomeryForm&) = delete;
    MontgomeryForm& operator=(const MontgomeryForm&) = delete;

    // Returns the largest valid modulus allowed for the constructor.
    static constexpr T max_modulus()
    {
        return (MontyType::max_modulus() >
            static_cast<U>(modular_arithmetic::ma_numeric_limits<T>::max()))
            ? ((modular_arithmetic::ma_numeric_limits<T>::max() == 0)
                ? static_cast<T>(modular_arithmetic::ma_numeric_limits<T>::max()
                                 - 1)
                : modular_arithmetic::ma_numeric_limits<T>::max())
            : static_cast<T>(MontyType::max_modulus());
    }

    // Returns the modulus given to the constructor
    T getModulus() const { return static_cast<T>(impl.getModulus()); }

    // Returns the converted value of the standard number 'a' into monty form.
    // Requires 0 <= a < modulus.
    MontgomeryValue convertIn(T a) const
    {
        HPBC_PRECONDITION(a >= 0);
        HPBC_PRECONDITION(a < static_cast<T>(impl.getModulus()));
        return impl.convertIn(static_cast<U>(a));
    }
    // Converts (montgomery value) x into a "normal" number; returns the result.
    // Guarantees 0 <= result < modulus.
    T convertOut(MontgomeryValue x) const
    {
        T a = static_cast<T>(impl.convertOut(x));
        HPBC_POSTCONDITION(a >= 0);
        HPBC_POSTCONDITION(a < static_cast<T>(impl.getModulus()));
        return a;
    }

    // Returns a unique (canonical) value representing the equivalence class of
    // x modulo the modulus.  You can not directly compare MontgomeryValues, but
    // you can call getCanonicalValue(), and then use standard equality or
    // inequality operators to compare the resulting CanonicalValues.
    CanonicalValue getCanonicalValue(MontgomeryValue x) const
    {
        MontgomeryValue ret = impl.getCanonicalValue(x);
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value 1.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(1)))
    // but it's more efficient (essentially zero cost) and more convenient.
    CanonicalValue getUnityValue() const
    {
        MontgomeryValue ret = impl.getUnityValue();
        HPBC_POSTCONDITION(impl.isCanonical(ret));
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value 0.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(0)))
    // but it's more efficient (essentially zero cost) and more convenient.
    CanonicalValue getZeroValue() const
    {
        MontgomeryValue ret = impl.getZeroValue();
        HPBC_POSTCONDITION(impl.isCanonical(ret));
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value
    // modulus-1 (which equals -1 (mod modulus)).  The call is equivalent to
    // getCanonicalValue(convertIn(static_cast<T>(modulus - 1))), but it's more
    // efficient (essentially zero cost) and more convenient.
    CanonicalValue getNegativeOneValue() const
    {
        MontgomeryValue ret = impl.getNegativeOneValue();
        HPBC_POSTCONDITION(impl.isCanonical(ret));
        return CanonicalValue(ret);
    }

    // Returns the modular product of (the montgomery values) x and y.
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.multiply(x, y);
    }

    // Returns the modular product of (the montgomery value) x multiplied by x.
    MontgomeryValue square(MontgomeryValue x) const
    {
        return impl.multiply(x, x);
    }

    // Calculates and returns the modular exponentiation of the montgomery value
    // 'base' to the power of (the type T variable) 'exponent'.
    MontgomeryValue pow(MontgomeryValue base, T exponent) const
    {
        HPBC_PRECONDITION(exponent >= 0);
        // This is a slightly optimized version of Algorithm 14.76, from
        // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
        // See also: hurchalla/modular_arithmetic/internal/impl_modular_pow.h
        MontgomeryValue result;
        if (exponent & static_cast<T>(1))
            result = base;
        else
            result = impl.getUnityValue();
        while (exponent > static_cast<T>(1))
        {
            exponent = static_cast<T>(exponent >> static_cast<T>(1));
            base = impl.multiply(base, base);
            if (exponent & static_cast<T>(1))
                result = impl.multiply(result, base);
        }
        return result;
    }

    // Returns the modular sum of (the montgomery values) x and y.
    MontgomeryValue add(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.add(x, y);
    }

    // Returns the modular difference of (the montgomery values) x and y.  More
    // precisely, x minus y.
    MontgomeryValue subtract(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.subtract(x, y);
    }
};


}} // end namespace

#endif
