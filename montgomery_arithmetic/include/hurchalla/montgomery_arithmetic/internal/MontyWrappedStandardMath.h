
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/internal/MontgomeryValue.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "hurchalla/montgomery_arithmetic/internal/compiler_macros.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


template <typename T>
class MontyWrappedStandardMath final {
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    using V = MontgomeryValue<T>;
    T modulus_;
public:
    using template_param_type = T;
    using montvalue_type = V;

    explicit MontyWrappedStandardMath(T modulus) : modulus_(modulus)
    {
        precondition2(modulus > 0);
    }
    MontyWrappedStandardMath(const MontyWrappedStandardMath&) = delete;
    MontyWrappedStandardMath& operator=(const MontyWrappedStandardMath&)=delete;

    // intended for use in postconditions/preconditions
    FORCE_INLINE bool isCanonical(V x) const
    {
        return (0 <= x.get() && x.get() < modulus_);
    }
    FORCE_INLINE T getModulus() const
    {
        return modulus_;
    }
//    FORCE_INLINE bool isValid(V x) const { return isCanonical(x); }

    FORCE_INLINE V convertIn(T a) const
    {
        precondition2(0 <= a && a < modulus_);
        return V(a);
    }
    FORCE_INLINE T convertOut(V x) const
    {
        precondition2(isCanonical(x));
        T ret = x.get();
        postcondition2(0 <= ret && ret < modulus_);
        return ret;
    }

    FORCE_INLINE V getCanonicalForm(V x) const
    {
        precondition2(isCanonical(x));
        return x;
    }

    FORCE_INLINE V getUnityValue() const        { return V(static_cast<T>(1)); }
    FORCE_INLINE V getZeroValue() const         { return V(static_cast<T>(0)); }
    FORCE_INLINE V getNegativeOneValue() const  { return V(modulus_ -
                                                           static_cast<T>(1)); }

    FORCE_INLINE V multiply(V x, V y) const
    {
        precondition2(isCanonical(x));
        precondition2(isCanonical(y));
        namespace ma = hurchalla::modular_arithmetic;
        T result = ma::modular_multiplication_prereduced_inputs(x.get(),
                                                             y.get(), modulus_);
        postcondition2(isCanonical(V(result)));
        return V(result);
    }
    FORCE_INLINE V square(V x) const
    {
        return multiply(x, x);
    }
    FORCE_INLINE V add(V x, V y) const
    {
        precondition2(isCanonical(x));
        precondition2(isCanonical(y));
        namespace ma = hurchalla::modular_arithmetic;
        T result =
          ma::modular_addition_prereduced_inputs(x.get(), y.get(), modulus_);
        postcondition2(isCanonical(V(result)));
        return V(result);
    }
    FORCE_INLINE V subtract(V x, V y) const
    {
        precondition2(isCanonical(x));
        precondition2(isCanonical(y));
        namespace ma = hurchalla::modular_arithmetic;
        T result =
          ma::modular_subtraction_prereduced_inputs(x.get(), y.get(), modulus_);
        postcondition2(isCanonical(V(result)));
        return V(result);
    }
};


}} // end namespace

#endif
