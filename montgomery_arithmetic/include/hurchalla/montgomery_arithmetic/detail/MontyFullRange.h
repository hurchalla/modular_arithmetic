// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// Let the theoretical constant R = 2^(ut_numeric_limits<T>::digits).
template <typename T>
class MontyFullRange : public MontyCommonBase<MontyFullRange, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    using BC = MontyCommonBase<::hurchalla::detail::MontyFullRange, T>;
    using BC::n_;
    using BC::inv_n_;
    using V = typename BC::MontgomeryValue;
public:
    using montvalue_type = V;
    using uint_type = T;
    using MontyTag = FullrangeTag;

    explicit MontyFullRange(T modulus) : BC(modulus) {}
    MontyFullRange(const MontyFullRange&) = delete;
    MontyFullRange& operator=(const MontyFullRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return (ut_numeric_limits<T>::max() % 2 == 0) ?
                static_cast<T>(ut_numeric_limits<T>::max() - 1) :
                ut_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE T getExtendedModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        return x;
    }

    // For PTAG (Performance Tag), see optimization_tag_structs.h
    template <bool smallModulus, class PTAG>
    HURCHALLA_FORCE_INLINE typename std::enable_if<smallModulus, V>::type
    famul(V x, V y, V z, bool& isZero) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        HPBC_PRECONDITION2(y.get() < n_);
        HPBC_PRECONDITION2(z.get() < n_);
        // Since smallModulus == true, we require modulus < R/2 (which allows
        // us to optimize this function).
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(n_ < Rdiv2);

        // x+y won't overflow:  x<n_<R/2 and y<n_<R/2, so  x+y < R/2 + R/2 == R.
        T sum = static_cast<T>(x.get() + y.get());
        // As a precondition, REDC requires  sum*z < n_*R.  This will always
        // be satisfied:  we know  sum = x+y < R  and  z < n_.  Therefore
        // sum*z < R*n_ == n_*R.
        // Though we would now like to call BC::multiply(sum, z, isZero), we
        // can't because sum<R does not satisfy its precondition of sum < n_.
        // Instead we mostly replicate BC::multiply() here, and we know from
        // above that we will be able to safely call REDC.
        T u_lo;
        T u_hi = unsigned_multiply_to_hilo_product(&u_lo, sum, z.get());
        // u_hi < n  implies that  sum*z == u < n*R.  See REDC_non_finalized()
        // in impl_REDC.h for proof.
        HPBC_ASSERT2(u_hi < n_);
        T result = REDC(u_hi, u_lo, n_, inv_n_, isZero, MontyTag(), PTAG());

        HPBC_POSTCONDITION2(isZero ==
              (getCanonicalValue(V(result)).get() == BC::getZeroValue().get()));
        HPBC_POSTCONDITION2(result < n_);
        return V(result);
    }
    template <bool smallModulus, class PTAG>
    HURCHALLA_FORCE_INLINE typename std::enable_if<!smallModulus, V>::type
    famul(V x, V y, V z, bool& isZero) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        HPBC_PRECONDITION2(y.get() < n_);
        HPBC_PRECONDITION2(z.get() < n_);
        // Since smallModulus==false, we have no extra requirements on modulus.
        // Unfortunately this means we can't do a simple non-modular add of
        // (x+y) prior to the multiply by z, because the sum would not
        // necessarily satisfy REDC's precondition that sum*z < R*n_.  REDC is
        // required for the montgomery multiply that this function must use.  To
        // see why we would be unable to guarantee the precondition:
        // We have  sum*z == (x+y)*z < (n_+n_)*n_ == 2*n_*n_.  Since n_ is type
        // T, we know n_ < R.  Thus we have  sum*z == 2*n_*n_ < 2*n_*R.  This
        // result of  sum*z < 2*n_*R  is insufficient to guarantee that we would
        // satisfy REDC's precondition requirement of  sum*z < n_*R.
        //
        // Instead we will use modular addition to get the sum, which results in
        // sum < n_ (thus guaranteeing  sum*z < n_*R).
        // In the end this version of famul(), which has the template argument
        // smallModulus==false, just wraps this class's add_canonical_value
        // and multiply functions.
        V sum = BC::add_canonical_value(x, y);
        V result = BC::multiply(sum, z, isZero, PTAG());

        HPBC_POSTCONDITION2(result.get() < n_);
        return result;
    }
};


}} // end namespace

#endif
