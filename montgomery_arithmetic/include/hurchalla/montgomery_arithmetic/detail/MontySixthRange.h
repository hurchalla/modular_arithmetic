// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SIXTH_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SIXTH_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/Redc.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


// Let the theoretical constant R = 2^(ut_numeric_limits<T>::digits).
template <typename T>
class MontySixthRange : public MontyCommonBase<MontySixthRange, T> {
    static_assert(util::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(util::ut_numeric_limits<T>::is_signed), "");
    static_assert(util::ut_numeric_limits<T>::is_modulo, "");
    using BC = MontyCommonBase<
                        ::hurchalla::montgomery_arithmetic::MontySixthRange, T>;
    using BC::n_;
    using BC::inv_n_;
    using V = typename BC::MontgomeryValue;
public:
    using montvalue_type = V;
    using template_param_type = T;
    using MontyTag = SixthrangeTag;

    explicit MontySixthRange(T modulus) : BC(modulus)
    {
        // MontySixthRange requires  modulus < R/6
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                      (util::ut_numeric_limits<T>::digits - 1));
        T Rdiv6 = static_cast<T>(Rdiv2 / 3);
        HPBC_PRECONDITION2(modulus < Rdiv6);
    }
    MontySixthRange(const MontySixthRange&) = delete;
    MontySixthRange& operator=(const MontySixthRange&) = delete;

    static constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                               (util::ut_numeric_limits<T>::digits - 1))/3 - 1);
    }

    HURCHALLA_FORCE_INLINE T getExtendedModulus() const
    {
        return static_cast<T>(2*n_);
    }

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        T cv = (x.get() < n_) ? x.get() : static_cast<T>(x.get() - n_);
        HPBC_POSTCONDITION2(cv < n_);
        return V(cv);
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V famul(V x, V y, V z, bool& isZero, PTAG) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < n_);   // y must be canonical
        HPBC_PRECONDITION2(z.get() < 2*n_);

        // x+y won't overflow:  we know  x < 2*n_ < 2*R/6  and  y < n_ < R/6,
        // so  x + y < 2*R/6 + R/6 == R/2.
        T sum = static_cast<T>(x.get() + y.get());
        // As a precondition, REDC requires  sum*z < n_*R.  This will always
        // be satisfied:  we know  sum = x+y < R/2  and  z < 2*n_.  Therefore
        // sum*z < (R/2)*(2*n_) == n_*R.
        // Though we would now like to call BC::multiply(sum, z, isZero), we
        // can't because sum<R/2 does not satisfy its precondition of sum < n_.
        // Instead we mostly replicate BC::multiply()) here, and we know from
        // above that we will be able to safely call REDC
        T u_lo;
        T u_hi = unsigned_multiply_to_hilo_product(&u_lo, sum, z.get());
        // u_hi < n  implies that  sum*z == u < n*R.  See REDC_non_finalized()
        // in Redc.h for proof.
        HPBC_ASSERT2(u_hi < n_);
        T result = REDC(u_hi, u_lo, n_, inv_n_, SixthrangeTag(), PTAG());
        isZero = isZeroRedcResult(result, n_, SixthrangeTag());

        HPBC_POSTCONDITION2(isZero ==
              (getCanonicalValue(V(result)).get() == BC::getZeroValue().get()));
        HPBC_POSTCONDITION2(result < 2*n_);  //REDC guarantees for SixthrangeTag
        return V(result);
    }
};


}} // end namespace

#endif
