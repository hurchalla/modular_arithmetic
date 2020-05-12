
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/monty_common.h"
#include "hurchalla/montgomery_arithmetic/internal/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "hurchalla/montgomery_arithmetic/internal/compiler_macros.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


template <typename T>
class MontyFullRange final : public MontyCommonBase<MontyFullRange, T> {
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    using MontyCommonBase<MontyFullRange, T>::n_;
    using MontyCommonBase<MontyFullRange, T>::neg_inv_n_;
    using typename MontyCommonBase<MontyFullRange, T>::V;
public:
    using montvalue_type = V;
    using template_param_type = T;

    explicit MontyFullRange(T modulus) :
                                  MontyCommonBase<MontyFullRange, T>(modulus) {}
    MontyFullRange(const MontyFullRange&) = delete;
    MontyFullRange& operator=(const MontyFullRange&) = delete;

    FORCE_INLINE bool isValid(V x) const { return (x.get() < n_); }
//    FORCE_INLINE bool isReduced(T a) const  { return (a < n_); }

    FORCE_INLINE T convertOut(V x) const
    {
        precondition2(x.get() < n_);
        T a = montout_non_minimized(x.get(), n_, neg_inv_n_);
        // montout_non_minimized() postconditions guarantee that since x < n_,
        // a < n_.  Thus 'a' is already minimized.
        T minimized_result = a;
        postcondition2(minimized_result < n_);
        return minimized_result;
    }

    FORCE_INLINE V getCanonicalForm(V x) const
    {
        precondition2(x.get() < n_);
        return x;
    }

    FORCE_INLINE V multiply(V x, V y) const
    {
        precondition2(x.get() < n_);
        precondition2(y.get() < n_);
        // x<n with y<n  will satisfy montmul_non_minimized's precondition
        // requirement that x*y < n*R.
        bool ovf;
        T prod = montmul_non_minimized(ovf, x.get(), y.get(), n_, neg_inv_n_);
        // montmul_non_minimized() postconditions guarantee the following
        T minimized_result = (ovf || prod >= n_) ? (prod - n_) : prod;
        postcondition2(minimized_result < n_);
        return V(minimized_result);
    }

    FORCE_INLINE V add(V x, V y) const
    {
        precondition2(x.get() < n_);
        precondition2(y.get() < n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_addition_prereduced_inputs(x.get(), y.get(), n_);
        postcondition2(z < n_);
        return V(z);
    }

    FORCE_INLINE V subtract(V x, V y) const
    {
        precondition2(x.get() < n_);
        precondition2(y.get() < n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_subtraction_prereduced_inputs(x.get(), y.get(), n_);
        postcondition2(z < n_);
        return V(z);
    }
};


}} // end namespace

#endif
