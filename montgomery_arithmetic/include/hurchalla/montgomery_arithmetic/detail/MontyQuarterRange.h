
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


// Let the theoretical constant R = 2^(ma_numeric_limits<T>::digits).
template <typename T>
class MontyQuarterRange final : public MontyCommonBase<MontyQuarterRange, T> {
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    using MontyCommonBase<MontyQuarterRange, T>::n_;
    using MontyCommonBase<MontyQuarterRange, T>::neg_inv_n_;
    using typename MontyCommonBase<MontyQuarterRange, T>::V;
public:
    using montvalue_type = V;
    using template_param_type = T;

    explicit MontyQuarterRange(T modulus) :
                                  MontyCommonBase<MontyQuarterRange, T>(modulus)
    {
        // MontyQuarterRange requires  modulus < R/4
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(modulus < Rdiv4);
    }

    MontyQuarterRange(const MontyQuarterRange&) = delete;
    MontyQuarterRange& operator=(const MontyQuarterRange&) = delete;

    static constexpr T max_modulus()
    {
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        return static_cast<T>(Rdiv4 - 1);
    }

    HURCHALLA_FORCE_INLINE bool isValid(V x) const { return (x.get() < 2*n_); }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        T a = montout_non_minimized(x.get(), n_, neg_inv_n_);
        // montout_non_minimized() postconditions guarantee the following
        T minimized_result = (a >= n_) ? static_cast<T>(a - n_) : a;
        HPBC_POSTCONDITION2(minimized_result < n_);
        return minimized_result;
    }

    HURCHALLA_FORCE_INLINE V getCanonicalForm(V x) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        T cf = (x.get() < n_) ? x.get() : static_cast<T>(x.get() - n_);
        HPBC_POSTCONDITION2(cf < n_);
        return V(cf);
    }

    HURCHALLA_FORCE_INLINE V multiply(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < 2*n_);
        // Since x<2*n and y<2*n, we know x*y < 4*n*n, and since we have a class
        // precondition that our modulus n < R/4, we know  x*y < 4*n*R/4 == n*R.
        // This satisfies montmul_non_minimized's precondition of x*y < n*R.
        bool ovf;
        T prod = montmul_non_minimized(ovf, x.get(), y.get(), n_, neg_inv_n_);

        // Since our constructor required modulus n < R/4, the postconditions of
        // montmul_non_minimized() guarantee  prod < 2*n.
        HPBC_POSTCONDITION2(prod < 2*n_);
        return V(prod);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < 2*n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_addition_prereduced_inputs(x.get(), y.get(), 
                                                          static_cast<T>(2*n_));
        HPBC_POSTCONDITION2(z < 2*n_);
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < 2*n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_subtraction_prereduced_inputs(x.get(), y.get(),
                                                          static_cast<T>(2*n_));
        HPBC_POSTCONDITION2(z < 2*n_);
        return V(z);
    }
};


}} // end namespace

#endif
