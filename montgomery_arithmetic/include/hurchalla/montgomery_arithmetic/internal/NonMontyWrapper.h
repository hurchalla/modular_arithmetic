
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_NON_MONTY_WRAPPER_H__INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_NON_MONTY_WRAPPER_H__INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/montycommonbase.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic {


template <typename T>
struct NonMontyWrapper final : public MontyCommonBase<NonMontyWrapper<T>> {
private:
    typedef MontgomeryValue<T> V;
    static_assert(std::is_unsigned<T>::value, "");
public:
    explicit NonMontyWrapper(T modulus) : modulus_(modulus) {}
    NonMontyWrapper(const NonMontyWrapper&) = delete;
    NonMontyWrapper& operator=(const NonMontyWrapper&) = delete;

    V convertIn(T a) const { return V(a); }
    T convertOut(V x) const
    {
        T ret = x.get();
        postcondition(ret < modulus_);
        return ret;
    }

    V getCanonicalForm(V x) const { return x; }

    V getUnityValue() const          { return convertIn(static_cast<T>(1)); }
    V getZeroValue() const           { return convertIn(static_cast<T>(0)); }
    V getNegativeOneValue() const    { return modulus_ - getUnityValue(); }

    V multiply(V x, V y) const
    {
        return V(modular_multiplication_prereduced_inputs(x.get(), y.get(), modulus_));
    }

    V add(V abar, V bbar) const noexcept      { return V(modular_addition_prereduced_inputs(abar.get(), bbar.get(), modulus_)); }
    V subtract(V abar, V bbar) const noexcept { return V(modular_subtraction_prereduced_inputs(abar.get(), bbar.get(), modulus_)); }



    bool isValid(V x) const { return (x.get() < modulus_); }
};
    // Returns the modular sum of (the montgomery values) x and y.  The return
    // value is in montgomery form but might not be canonical - call
    // getCanonicalForm() to use it in comparisons.
    V add(V x, V y) const noexcept { return impl.add(x, y); }

    // Returns the modular difference of (the montgomery values) x and y.  More
    // precisely, x minus y.  The return value is in montgomery form but might
    // not be canonical - call getCanonicalForm() to use it in comparisons.
    V subtract(V x, V y) const noexcept { return impl.subtract(x, y); }


}} // end namespace

#endif
