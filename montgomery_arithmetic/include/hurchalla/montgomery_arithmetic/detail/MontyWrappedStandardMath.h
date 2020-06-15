
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/detail/MontgomeryValue.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"

namespace hurchalla { namespace montgomery_arithmetic {


// This class provides a standard modular arithmetic implementation, wrapped
// inside a Monty template.  This allows standard modular arithmetic to be used
// in a generic MontgomeryForm instantation.
template <typename T>
class MontyWrappedStandardMath final {
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    using V = MontgomeryValue<T>;
    T modulus_;
public:
    using template_param_type = T;
    using montvalue_type = V;

    explicit MontyWrappedStandardMath(T modulus) : modulus_(modulus)
    {
        HPBC_PRECONDITION2(modulus > 0);
    }
    MontyWrappedStandardMath(const MontyWrappedStandardMath&) = delete;
    MontyWrappedStandardMath& operator=(const MontyWrappedStandardMath&)=delete;

    static constexpr T max_modulus()
    {
        return modular_arithmetic::ma_numeric_limits<T>::max();
    }

    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        return (0 <= x.get() && x.get() < modulus_);
    }
    HURCHALLA_FORCE_INLINE T getModulus() const
    {
        return modulus_;
    }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        HPBC_PRECONDITION2(0 <= a && a < modulus_);
        return V(a);
    }
    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        T ret = x.get();
        HPBC_POSTCONDITION2(0 <= ret && ret < modulus_);
        return ret;
    }

    HURCHALLA_FORCE_INLINE V getCanonicalForm(V x) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        return x;
    }

    HURCHALLA_FORCE_INLINE V getUnityValue() const
    {
        return V(static_cast<T>(1));
    }
    HURCHALLA_FORCE_INLINE V getZeroValue() const
    {
        return V(static_cast<T>(0));
    }
    HURCHALLA_FORCE_INLINE V getNegativeOneValue() const 
    {
        return V(static_cast<T>(modulus_ - static_cast<T>(1)));
    }

    HURCHALLA_FORCE_INLINE V multiply(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        namespace ma = hurchalla::modular_arithmetic;
        T result = ma::modular_multiplication_prereduced_inputs(x.get(),
                                                             y.get(), modulus_);
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        namespace ma = hurchalla::modular_arithmetic;
        T result =
          ma::modular_addition_prereduced_inputs(x.get(), y.get(), modulus_);
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        namespace ma = hurchalla::modular_arithmetic;
        T result =
          ma::modular_subtraction_prereduced_inputs(x.get(), y.get(), modulus_);
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
};


}} // end namespace

#endif
