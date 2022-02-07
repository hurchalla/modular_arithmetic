// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// The name "Quarterrange" signifies that the modulus must be less than R/4,
// where  R = 2^(ut_numeric_limits<T>::digits).  For example, if T is uint64_t
// then R = 2^64 and R/4 == 2^62, and thus it would require  modulus < 2^62.

// The MontyQuarterRange class functions require/allow an unusual input range:
// for an input x, they allow  0 <= x < 2*n, where n is the modulus.  Similarly,
// the return value range will be  0 <= returnValue < 2*n.  Obviously neither
// inputs nor outputs necessarily belong to the minimal residue class modulo n -
// i.e. they might not be fully reduced, modulo n.  Note that the algorithm for
// montgomery REDC requires that  u = x*y < n*R;  this will always be satisfied
// for any multiplication x*y of Quarterrange montgomery values.  To see why,
// keep in mind that Quarterrange requires n < R/4  and that all inputs are less
// than 2*n.  Thus the multiplication
// u = x*y < (2*n)*(2*n) == (4*n)*n < (4*n)*(R/4) == n*R,  which means u < n*R,
// as required.  For more details, see also section 5 from
// "Montgomery's Multiplication Technique: How to Make It Smaller and Faster"
// https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps


struct TagMontyQuarterrange final {};  // IDs MontyQuarterRange independent of T


// struct used internally by MontyQuarterRange
template <typename T>
struct MontyQRValueTypes {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    // regular montgomery value type
    struct V : public BaseMontgomeryValue<T> {
        HURCHALLA_FORCE_INLINE V() = default;
     protected:
        template <typename> friend class MontyQuarterRange;
        HURCHALLA_FORCE_INLINE explicit V(T a) : BaseMontgomeryValue<T>(a) {}
    };
    // canonical montgomery value type
    struct C : public V {
        HURCHALLA_FORCE_INLINE C() = default;
        HURCHALLA_FORCE_INLINE friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        HURCHALLA_FORCE_INLINE friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }
     protected:
        template <template<class> class, template<class> class, typename>
          friend class MontyCommonBase;
        template <typename> friend class MontyQuarterRange;
        HURCHALLA_FORCE_INLINE explicit C(T a) : V(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    using FV = C;
};


// Let the theoretical constant R = 2^(ut_numeric_limits<T>::digits).
template <typename T>
class MontyQuarterRange final : public
                      MontyCommonBase<MontyQuarterRange, MontyQRValueTypes, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    static_assert(ut_numeric_limits<T>::digits >= 2, "");
    using BC = MontyCommonBase<::hurchalla::detail::MontyQuarterRange,
                               ::hurchalla::detail::MontyQRValueTypes, T>;
    using BC::n_;
    using typename BC::V;
    using typename BC::C;
    using FV = typename MontyQRValueTypes<T>::FV;
 public:
    using MontyTag = TagMontyQuarterrange;
    using uint_type = T;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;

    using BC::fmadd;
    using BC::fmsub;

    explicit MontyQuarterRange(T modulus) : BC(modulus)
    {
        // MontyQuarterRange requires  modulus < R/4
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(modulus < Rdiv4);
    }
    MontyQuarterRange(const MontyQuarterRange&) = delete;
    MontyQuarterRange& operator=(const MontyQuarterRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 2)) - 1);
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        // this static_assert guarantees 0 <= x.get()
        static_assert(!(ut_numeric_limits<decltype(x.get())>::is_signed), "");
        HPBC_PRECONDITION2(x.get() < static_cast<T>(2*n_));
        T c = static_cast<T>(x.get() - n_);
        HURCHALLA_CMOV(x.get() < n_, c, x.get()); //if x.get()<n_, set c=x.get()
        HPBC_POSTCONDITION2(c < n_);
        return C(c);
    }

    static_assert(std::is_same<FV, C>::value, "");
    // Note: fmsub and fmadd with FusingValue (FV) arguments will match to
    // fmsub and fmadd with CanonicalValue args, since C is_same as FV.

    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        return getCanonicalValue(x);
    }

    using BC::add;
    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T n2 = static_cast<T>(2*n_);
        T result = modular_addition_prereduced_inputs(x.get(), y.get(), n2);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V add(V x, C y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(y.get() < n_);
        T result = mont_add_canonical_value<T>::call(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    using BC::subtract;
    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T n2 = static_cast<T>(2*n_);
        T result = modular_subtraction_prereduced_inputs(x.get(), y.get(), n2);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract(V x, C y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(y.get() < n_);
        T result = mont_subtract_canonical_value<T>::call(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: subtract(C, V) will match to subtract(V x, V y) above

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T result = absolute_value_difference(x.get(), y.get());
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: unordered_subtract(C, V) and unordered_subtract(V, C) will match to
    // unordered_subtract(V x, V y) above

private:
    // functions called by the 'curiously recurring template pattern' base (BC).
    friend BC;

    // note: PTAG gets ignored here
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, PTAG) const
    {
        HPBC_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        bool isNegative;  // ignored
        T result = REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
        resultIsZero = (result == 0);
        T sum = static_cast<T>(result + n_);
        HPBC_POSTCONDITION2(0 < sum && sum < static_cast<T>(2*n_));
        return V(sum);
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(T u_hi, T u_lo, PTAG) const
    {
        bool resultIsZero;  // ignored
        return montyREDC(resultIsZero, u_hi, u_lo, PTAG());
    }
    // return the high word of the product, and write the low word of the
    // product to u_lo.
    HURCHALLA_FORCE_INLINE T multiplyToHiLo(T& u_lo, V x, V y) const
    {
        return unsigned_multiply_to_hilo_product(u_lo, x.get(), y.get());
    }
    HURCHALLA_FORCE_INLINE T squareToHiLo(T& u_lo, V x) const
    {
        return unsigned_multiply_to_hilo_product(u_lo, x.get(), x.get());
    }
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        return (x.get() < static_cast<T>(2*n_));
    }
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        return (0 <= x.get() && x.get() < n_);
    }
    // get a Natural number (i.e. number >= 0) congruent to x (mod n)
    HURCHALLA_FORCE_INLINE T getNaturalEquivalence(V x) const
    {
        return x.get();
    }
};


}} // end namespace

#endif
