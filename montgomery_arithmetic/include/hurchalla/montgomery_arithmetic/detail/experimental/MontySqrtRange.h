// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SQRT_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SQRT_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/experimental/platform_specific/montadd_sqrt_range.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/platform_specific/montsub_sqrt_range.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/negative_inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/detail/BaseMontgomeryValue.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#  pragma warning(disable : 4309)
#endif

namespace hurchalla { namespace detail {


// For discussion purposes, let R = 1<<(ut_numeric_limits<T>::digits).  For
// example if T is uint64_t, then R = 1<<64.

// MontySqrtRange uses optimizations based on input and output V values being
// 0 < val <= n_, and on modulus < sqrtR.
// These restrictions allow us to implement a more efficient version of the REDC
// algorithm  (in the function msr_REDC_non_minimized()), by omitting some
// conditionals and calculations that would normally be needed.


struct TagMontySqrtrange final {};  //identifies MontySqrtRange independent of T


// The class member variable names are based on the webpage
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
template <typename T>
class MontySqrtRange final {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    const T n_;   // the modulus
    const T r_mod_n_;
    const T neg_inv_n_;
    const T r_squared_mod_n_;

    struct V : public BaseMontgomeryValue<T> {  // regular montgomery value type
        HURCHALLA_FORCE_INLINE V() = default;
     protected:
        friend MontySqrtRange;
        HURCHALLA_FORCE_INLINE explicit V(T a) : BaseMontgomeryValue<T>(a) {}
    };
    struct C : public V {                     // canonical montgomery value type
        HURCHALLA_FORCE_INLINE C() = default;
        HURCHALLA_FORCE_INLINE friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        HURCHALLA_FORCE_INLINE friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }
     protected:
        friend MontySqrtRange;
        HURCHALLA_FORCE_INLINE explicit C(T a) : V(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    struct FV : public V {
        HURCHALLA_FORCE_INLINE FV() = default;
     protected:
        friend MontySqrtRange;
        HURCHALLA_FORCE_INLINE explicit FV(T a) : V(a) {}
    };

 public:
    using MontyTag = TagMontySqrtrange;

    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;
    using uint_type = T;

    explicit MontySqrtRange(T modulus) :
                n_(modulus),
                r_mod_n_(getRModN(n_)),
                neg_inv_n_(negative_inverse_mod_R(n_)),
                r_squared_mod_n_(modular_multiplication_prereduced_inputs(
                                                        r_mod_n_, r_mod_n_, n_))
    {
        static constexpr int bitsT = ut_numeric_limits<T>::digits;
        static_assert(bitsT % 2 == 0, "");   // bitsT divisible by 2
        // MontySqrtRange requires  modulus < sqrt(R)
        static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);
        HPBC_PRECONDITION2(1 < modulus && modulus < sqrtR);
        HPBC_PRECONDITION2(modulus % 2 == 1);

        // Note: unityValue == (the montgomery form of 1)==(1*R)%n_ == r_mod_n_.
        // getRModN() guarantees the below.  getUnityValue() and
        // getNegativeOneValue() rely on it.
        HPBC_INVARIANT2(0 < r_mod_n_ && r_mod_n_ < n_);
        // Since n_ == modulus is odd and n_ > 1, n_ can not divide R*R because
        // R is a power of 2.
        // Thus  r_squared_mod_n_ == R*R (mod n_) != 0.  convertIn relies on it.
        HPBC_INVARIANT2(0 < r_squared_mod_n_ && r_squared_mod_n_ < n_);
    }
    MontySqrtRange(const MontySqrtRange&) = delete;
    MontySqrtRange& operator=(const MontySqrtRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        static_assert(ut_numeric_limits<T>::digits%2 == 0, "");
        return (static_cast<T>(1) << (ut_numeric_limits<T>::digits/2)) - 1;
    }

 private:
    static T getRModN(T n)
    {
        HPBC_PRECONDITION2(n % 2 == 1);
        HPBC_PRECONDITION2(n > 1);
        // Assign a tmp T variable rather than directly using the intermediate
        // expression, in order to avoid a negative value (and a wrong answer)
        // in cases where 'n' would be promoted to type 'int'.
        T tmp = static_cast<T>(static_cast<T>(0) - n);
        // Compute R%n.  For example if R==1<<64, arithmetic wraparound behavior
        // of the unsigned integral type T results in (0 - n) representing
        // ((1<<64) - n).  Thus,
        // rModN = R%n == (1<<64)%n == ((1<<64) - n)%n == (0-n)%n
        T rModN = static_cast<T>(tmp % n);
        // Since n is odd and > 1, n does not divide R because R is a power of
        // 2.  Thus, rModN != 0
        HPBC_POSTCONDITION2(0 < rModN && rModN < n);
        return rModN;
    }

    // This function is based on REDC_non_minimized() in https://github.com/hurchalla/modular_arithmetic/blob/66281af1639031b04bdaf9b916e5d5638d3ded25/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/platform_specific/RedcLargeR.h#L44
    // It is an adaptation of REDC_non_minimized that omits calculations that
    // are not needed, given this function's preconditions of n < sqrt(R), and
    // u < R (i.e. u_hi == 0).  The precondition of u_hi == 0 is expressed
    // simply by the lack of a u_hi parameter; u_hi is implicitly treated as
    // zero inside this function.
    HURCHALLA_FORCE_INLINE T msr_REDC_non_minimized(T u_lo) const
    {
        // For casts, we want to use types that are protected from surprises and
        // undefined behavior caused by unsigned integral promotion rules in C++
        // https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/
        using P = typename safely_promote_unsigned<T>::type;
        static_assert(ut_numeric_limits<P>::is_modulo, "");
        HPBC_PRECONDITION2(u_lo != 0);
        // Implicitly, u_hi == 0.  And thus u = (u_hi*R + u_lo) == u_lo < R.
        // Since this class guarantees n_ > 1, we know u < R < n_*R, which
        // satisfies the basic requirement of montgomery REDC that u < n_*R.

        // assert(n_ * neg_inv_n_ ≡ -1 (mod R))
        HPBC_PRECONDITION2(
              static_cast<T>(static_cast<P>(n_) * static_cast<P>(neg_inv_n_)) ==
              static_cast<T>(static_cast<P>(0) - static_cast<P>(1))
              );

        // compute  m = (u * neg_inv_n_) % R
        T m = static_cast<T>(static_cast<P>(u_lo) * static_cast<P>(neg_inv_n_));
        T mn_lo;
        T mn_hi = unsigned_multiply_to_hilo_product(mn_lo, m, n_);
        // mn = m*n_.  Since m=(u_lo*neg_inv_n_)%R, we know m < R, and thus
        // mn < R*n_.  Therefore mn == mn_hi*R + mn_lo < R*n_, and
        // mn_hi*R < R*n_ - mn_lo <= R*n_, and thus  mn_hi < n_.
        // *** Assertion #1 ***
        HPBC_ASSERT2(mn_hi < n_);
        // compute t_hi = (u_hi + mn_hi) % R.  Since we know u_hi == 0, we
        // simply omit the addition of u_hi.
        T t_hi = mn_hi;
        // The REDC algorithm guarantees (u_lo + mn_lo) % R == 0.
        HPBC_ASSERT2(static_cast<T>(u_lo + mn_lo) == static_cast<T>(0));
        // REDC_non_minimized() would normally next calculate
        // t_hi += (u_lo != 0);
        // However, we have know by precondition that u_lo != 0.  And so the
        // calculation of t_hi simplifies to
        t_hi = static_cast<T>(t_hi + static_cast<T>(1));
        // REDC_non_minimized() would normally next calculate
        // ovf = (t_hi < u_hi);
        // But we know u_hi == 0, so ovf = (t_hi < u_hi) == (t_hi < 0) == false.
        // Thus  ovf = false.

        // The discussion prior to Assertion #1 proves that mn_hi < n_.
        // Therefore, 0 < mn_hi + 1 < n_ + 1.  Since t_hi = mn_hi + 1, we know
        HPBC_POSTCONDITION2(0 < t_hi && t_hi <= n_);
        // From REDC_non_minimized()'s Postcondition #1, we know
        //   T minimized_result = (ovf || t_hi >= n_) ? (t_hi - n_) : t_hi;
        //   HPBC_POSTCONDITION2(minimized_result < n_);
        // Since  ovf == false  and  0 < t_hi <= n_,  we can simplify this to
        if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
            T minimized_result = (t_hi == n_) ? static_cast<T>(0) : t_hi;
            HPBC_POSTCONDITION2(minimized_result < n_);
        }
        // From REDC_non_minimized() Postcondition #3, we know
        //    HPBC_POSTCONDITION2((u_hi == 0 && u_lo < n_) ? t_hi < n_ : true);
        // Since u_hi == 0, we can simplify this to
        HPBC_POSTCONDITION2((u_lo < n_) ? t_hi < n_ : true);
        // return the non-minimized result
        return t_hi;
    }

    HURCHALLA_FORCE_INLINE T msr_montmul_non_minimized(T x, T y) const
    {
        // As in msr_REDC_non_minimized(), protect against undefined behavior:
        using P = typename safely_promote_unsigned<T>::type;
        static_assert(ut_numeric_limits<P>::is_modulo, "");
        static constexpr int bit_width_T = ut_numeric_limits<T>::digits;
        static_assert(bit_width_T % 2 == 0, "");   // bit_width_T divisible by 2
        static constexpr T sqrtR = static_cast<T>(1) << (bit_width_T / 2);
        HPBC_PRECONDITION2(0 < x && x < sqrtR);
        HPBC_PRECONDITION2(0 < y && y < sqrtR);

        // Since x < sqrtR and y < sqrtR,  x*y < sqrtR*sqrtR == R.
        // We have x*y < R, so x*y will fit in type T without overflow.
        T u_lo = static_cast<T>(static_cast<P>(x) * static_cast<P>(y));
        T result = msr_REDC_non_minimized(u_lo);

        HPBC_POSTCONDITION2(0 < result && result <= n_);
        return result;
    }

    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        return (0 < x.get() && x.get() <= n_);
    }

    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        C cfx = getCanonicalValue(x);
        bool good = isValid(x);
        return (x.get() == cfx.get() && good);
    }

 public:
    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    // We require a < sqrtR, which is a bit of a hack since MontgomeryForm class
    // expects that any T value >= 0 is ok to use as input for convertIn.
    // Ideally we would address this by renaming this class to something like
    // MontyDoubleWidth, and allowing all MontgomeryValues for this class to be
    // any T value >= 0, while setting
    // T2 = sized_uint<2* ut_numeric_limits<T>::digits>::type, and producing a
    // compile time error if
    // ut_numeric_limits<T2>::digits > HURCHALLA_TARGET_BIT_WIDTH.
    // n_, r_mod_n_, neg_inv_n_, r_squared_mod_n_  would all be type T2, and
    // most function calls in this class would work with type T2 values, and V
    // would wrap type T2.  While this all sounds like a big change, it is not-
    // this class effectively already works like this proposal, with the current
    // type T playing the role of the proposed T2, and convertIn's precondition
    // a < sqrtR playing a psuedo-role of the proposed T.
    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        static constexpr int bitsT = ut_numeric_limits<T>::digits;
        static_assert(bitsT % 2 == 0, "");
        static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);
        HPBC_PRECONDITION2(0 <= a && a < sqrtR);

        HPBC_INVARIANT2(1 < n_ && n_ < sqrtR);
        HPBC_INVARIANT2(0 < r_squared_mod_n_ && r_squared_mod_n_ < n_);
        // thus:  0 < r_squared_mod_n_ < sqrtR
        T result;
        if HURCHALLA_LIKELY(a > 0) {
            // We have  0 < a < sqrtR  and  0 < r_squared_mod_n_ < sqrtR,  which
            // satisfies the preconditions for msr_montmul_non_minimized().
            result = msr_montmul_non_minimized(a, r_squared_mod_n_);
        } else {
            HPBC_ASSERT2(a == 0);
            // We can't use msr_montmul_non_minimized() here, because it
            // requires nonzero inputs.  We must treat a == 0 as a special case:
            // a*R (mod n) ≡ 0*R (mod n) ≡ 0 (mod n) ≡ n (mod n).
            result = n_;
        }

        // both clauses of the if/else generate
        HPBC_POSTCONDITION2(0 < result && result <= n_);
        // Since 0 < result <= n, we don't want to reduce mod n;  result is in
        // the canonical form required by most of the class functions.
        return V(result);
    }

    HURCHALLA_FORCE_INLINE C getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_,
        // and 0 < r_mod_n_ < n_.
        HPBC_INVARIANT2(isCanonical(V(r_mod_n_)));
        return C(r_mod_n_);
    }

    HURCHALLA_FORCE_INLINE C getZeroValue() const
    {
        // We want returnVal == (0*R)%n_, but since isValid() requires
        // 0 < returnVal <= n_, we return n_ (n_ ≡ 0 (mod n_))
        HPBC_INVARIANT2(isCanonical(V(n_)));
        return C(n_);
    } 

    HURCHALLA_FORCE_INLINE C getNegativeOneValue() const
    {
        // We want to get returnVal = getCanonicalValue(subtract(getZeroValue(),
        //                                               getUnityValue())).
        //   getZeroValue() returns n_, and getUnityValue() returns r_mod_n_.
        //   Therefore the subtraction results in the equivalence class
        //   (n_ - r_mod_n_) (mod n_). The constructor established the invariant
        //   0 < r_mod_n_ < n_.  Thus we know  0 < n_ - r_mod_n_ < n_.  This
        //   means (n_ - r_mod_n_)  satisfies isValid() and getCanonicalValue().
        HPBC_INVARIANT2(n_ > r_mod_n_);
        T negOne = static_cast<T>(n_ - r_mod_n_);
        HPBC_ASSERT2(0 < negOne && negOne < n_);

        HPBC_INVARIANT2(isCanonical(V(negOne)));
        return C(negOne);
    }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);

        T prod = msr_REDC_non_minimized(x.get());

        // msr_REDC_non_minimized() postconditions guarantee the following
        HPBC_POSTCONDITION2(0 < prod && prod <= n_);

        T minimized_result;
        if HURCHALLA_LIKELY(prod != n_)
            minimized_result = prod;
        else
            minimized_result = static_cast<T>(0);
        HPBC_POSTCONDITION2(minimized_result < n_);
        return minimized_result;
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        return C(x.get());
    }

    HURCHALLA_FORCE_INLINE V negate(V x) const
    {
        return subtract(getZeroValue(), x);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        T a = x.get();
        T b = y.get();
        HPBC_PRECONDITION2(0 < a && a <= n_);
        HPBC_PRECONDITION2(0 < b && b <= n_);
        HPBC_INVARIANT2(n_ > 0);

        T result = montadd_sqrt_range<T>::call(a, b, n_);

        HPBC_POSTCONDITION2(0 < result && result <= n_);
        return V(result);
    }
    // Note: add(V, C) and add(C, V) will match to add(V x, V y) above.
    HURCHALLA_FORCE_INLINE C add(C x, C y) const
    {
        V v = add(V(x), V(y));
        // add(V,V) returns a result such that 0 < result <= n.  Thus our
        // result v satisfies:
        HPBC_POSTCONDITION2(isCanonical(v));
        return C(v.get());
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        T a = x.get();
        T b = y.get();
        HPBC_PRECONDITION2(0 < a && a <= n_);
        HPBC_PRECONDITION2(0 < b && b <= n_);
        HPBC_INVARIANT2(n_ > 0);

        T result = montsub_sqrt_range<T>::call(a, b, n_);

        HPBC_POSTCONDITION2(0 < result && result <= n_);
        return V(result);
    }
    // Note: subtract(V, C) and subtract(C, V) will match to subtract(V x, V y)
    // above.
    HURCHALLA_FORCE_INLINE C subtract(C x, C y) const
    {
        V v = subtract(V(x), V(y));
        // subtract(V,V) returns a result such that 0 < result <= n.  Thus our
        // result v satisfies:
        HPBC_POSTCONDITION2(isCanonical(v));
        return C(v.get());
    }

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        // we can't improve efficiency much over plain subtract,
        // so just delegate to subtract
        return subtract(x, y);
    }
    // Note: unordered_subtract(V, C) and unordered_subtract(C, V) will match
    // to unordered_subtract(V x, V y) above.


    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, bool& isZero, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        HPBC_PRECONDITION2(0 < y.get() && y.get() <= n_);

        // Since we know n < sqrtR  (guaranteed by the constructor),  and x < n
        // and  y < n,  we have  x < sqrtR  and  y < sqrtR,  which satisfies the
        // preconditions of msr_montmul_non_minimized().
        T prod = msr_montmul_non_minimized(x.get(), y.get());
        isZero = (getCanonicalValue(V(prod)).get() == getZeroValue().get());

        // msr_montmul_non_minimized() postconditions guarantee the following
        HPBC_POSTCONDITION2(0 < prod && prod <= n_);
        // Since 0 < prod <= n, we don't want to reduce mod n;  prod is in the
        // canonical form required by most of the class functions.
        return V(prod);
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, C z, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        HPBC_PRECONDITION2(0 < y.get() && y.get() <= n_);

        // Unfortunately for MontySqrtRange, it's not possible to get anything
        // more than perhaps a small efficiency advantage from a fused
        // multiply/add - in principle the small advantage could come from
        // inserting into (a copy of) msr_REDC_non_minimized() a modular add of
        // 'z' with 1 which occurs during the multiplications - thus the modular
        // add by 1 would not increase the latency.  The addition by 1 at the
        // end of (the copy of) msr_REDC_non_minimized() would be removed and
        // replaced with this->add(z_plus_one, REDC_result).  The REDC_result
        // would be 0 <= REDC_result < n, which is invalid for add(), so some
        // details would need to be worked out.
        // In the end this would decrease latency by 1 cycle compared to using
        // a combination of multiply() followed by add().  It would likely
        // increase the number of uops, which is not ideal.
        // Certainly for now, we will just implement fmadd as a wrapper of
        // multiply() followed by add().  If MontySqrtRange seems to be
        // beneficial enough, we can consider implementing the optimization
        // discussed here.

        bool isZero;
        V prod = multiply(x, y, isZero, PTAG());
        V sum = add(prod, z);

        HPBC_POSTCONDITION2(0 < sum.get() && sum.get() <= n_);
        // Since 0 < sum <= n, we don't want to reduce mod n;  sum is in the
        // canonical form required by most of the class functions.
        return sum;
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, C z, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        HPBC_PRECONDITION2(0 < y.get() && y.get() <= n_);

        // See the optimization discussion inside fmadd() - it applies here too.
        bool isZero;
        V prod = multiply(x, y, isZero, PTAG());
        V diff = subtract(prod, z);

        HPBC_POSTCONDITION2(0 < diff.get() && diff.get() <= n_);
        // Since 0 < diff <= n, we don't want to reduce mod n;  diff is in the
        // canonical form required by most of the class functions.
        return diff;
    }

    // Note: internal to MontySqrtRange, the contents of FusingValue (FV) and
    // CanonicalValue (C) variables are interchangeable.  Other Monty types
    // use FV and C as completely distinct types, and so for genericity we
    // always present C and FV to the outside world as being unrelated.
    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        C cv = getCanonicalValue(x);
        return FV(cv.get());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmadd(x, y, cv, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmsub(x, y, cv, PTAG());
    }

    // Returns the greatest common divisor of the standard representations
    // (non-montgomery) of both x and the modulus, using the supplied functor.
    // The functor must take two integral arguments of the same type and return
    // the gcd of those two arguments.  Usually you would make the functor's
    // operator() a templated function, where the template parameter is the
    // unknown type of the integral arguments.  Or more simply, you can just use
    // a lambda, with 'auto' type for the function parameters.
    template <class F>
    HURCHALLA_FORCE_INLINE T gcd_with_modulus(V x, const F& gcd_functor) const
    {
        HPBC_INVARIANT2(n_ > 0);
        // See the member function gcd_with_modulus() in MontyCommonBase.h for
        // proof that gcd(x.get(), n_) == gcd(convertOut(x), n_).
        // We want to return the value  q = gcd(convertOut(x), n_).
        // By relying on the equivalence of those two gcds, we can instead
        // compute and return p = gcd(x.get(), n_)  which we can compute more
        // efficiently than q.
        T p = gcd_functor(x.get(), n_);
        // Our postconditions assume the Functor implementation is correct.
        HPBC_POSTCONDITION2(0 < p && p <= n_ && (x.get() == 0 || p <= x.get()));
        HPBC_POSTCONDITION2(n_ % p == 0);
        HPBC_POSTCONDITION2(x.get() % p == 0);
        return p;
    }

    // This class doesn't do anything special for square functions.
    // It just delegates to the functions above.
    // --------
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V square(V x, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        bool isZero;
        return multiply(x, x, isZero, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareSub(V x, C cv, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        return fmsub(x, x, cv, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareAdd(V x, C cv, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        return fmadd(x, x, cv, PTAG());
    }
};


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
