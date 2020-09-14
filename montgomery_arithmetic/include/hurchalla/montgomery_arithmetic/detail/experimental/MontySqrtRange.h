// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SQRT_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SQRT_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/platform_specific/montadd_sqrt_range.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/platform_specific/montsub_sqrt_range.h"
#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#  pragma warning(disable : 4309)
#endif

namespace hurchalla { namespace montgomery_arithmetic {


// For discussion purposes, let R = 2^(ma_numeric_limits<T>::digits).  For
// example if T is uint64_t, then R = 2^64.
//
// This function is based on REDC_non_minimized() from RedcLargeR.h.  It is
// altered to omit calculations that are not needed, given its preconditions of
// n < sqrt(R), and u < R (i.e. u_hi == 0).  The precondition of u_hi == 0 is
// expressed simply by the lack of a u_hi parameter; u_hi is implicitly treated
// as zero inside this function.

template <typename T>
HURCHALLA_FORCE_INLINE T msr_REDC_non_minimized(T u_lo, T n, T neg_inv_n)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T>::is_modulo, "");

    // For casts, we want to use types that are protected from surprises and
    // undefined behavior due to the unsigned integral promotion rules in C++.
    // https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/
    using V = typename safely_promote_unsigned<T>::type;
    static_assert(ma::ma_numeric_limits<V>::is_modulo, "");

    HPBC_PRECONDITION2(n > 1);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(u_lo != 0);
    // Implicitly, u_hi == 0.  And thus u = (u_hi*R + u_lo) == u_lo < R.  Since
    // we have the precondition n > 1, u < R < n*R, which satisfies the basic
    // requirement of montgomery REDC that u < n*R.

    // assert(n * neg_inv_n ≡ -1 (mod R))
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<V>(n) * static_cast<V>(neg_inv_n)) ==
                static_cast<T>(static_cast<V>(0) - static_cast<V>(1))
                );

    // compute  m = (u * neg_inv_n) % R
    T m = static_cast<T>(static_cast<V>(u_lo) * static_cast<V>(neg_inv_n));

    T mn_lo;
    T mn_hi = unsigned_multiply_to_hilo_product(&mn_lo, m, n);

    // mn = m*n.  Since m=(u_lo*neg_inv_n)%R, we know m < R, and thus  mn < R*n.
    // Therefore mn == mn_hi*R + mn_lo < R*n, and mn_hi*R < R*n - mn_lo <= R*n,
    // and thus  mn_hi < n.
    // *** Assertion #1 ***
    HPBC_ASSERT2(mn_hi < n);

    // compute t_hi = (u_hi + mn_hi) % R.  Since we know u_hi == 0, we simply
    // omit the addition of u_hi.
    T t_hi = mn_hi;

    // The REDC algorithm guarantees (u_lo + mn_lo) % R == 0.
    HPBC_ASSERT2(static_cast<T>(u_lo + mn_lo) == static_cast<T>(0));
    // REDC_non_minimized() would normally next calculate
    // t_hi += (u_lo != 0);
    // However, we have know by precondition that u_lo != 0.  The calculation of
    // t_hi simplifies to
    t_hi = static_cast<T>(t_hi + static_cast<T>(1));

    // REDC_non_minimized() would normally next calculate
    // ovf = (t_hi < u_hi);
    // But we know u_hi == 0, so ovf = (t_hi < u_hi) == (t_hi < 0) == false.
    // Thus  ovf = false.

    // The discussion prior to Assertion #1 proves that mn_hi < n, and therefore
    // 0 < mn_hi + 1 < n + 1.  Since t_hi = mn_hi + 1, we know  0 < t_hi <= n.
    HPBC_POSTCONDITION2(0 < t_hi && t_hi <= n);
    // From REDC_non_minimized()'s Postcondition #1, we know
    //   T minimized_result = (ovf || t_hi >= n) ? (t_hi - n) : t_hi;
    //   HPBC_POSTCONDITION2(minimized_result < n);
    // Since  ovf == false  and  0 < t_hi <= n,  we can simplify this to
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (t_hi == n) ? static_cast<T>(0) : t_hi;
        HPBC_POSTCONDITION2(minimized_result < n);
    }

    // From REDC_non_minimized() Postcondition #3, we know
    //    HPBC_POSTCONDITION2((u_hi == 0 && u_lo < n) ? t_hi < n : true);
    // Since u_hi == 0, we can simplify this to
    HPBC_POSTCONDITION2((u_lo < n) ? t_hi < n : true);

    // return the non-minimized result
    return t_hi;
}


template <typename T>
HURCHALLA_FORCE_INLINE T msr_montmul_non_minimized(T x, T y, T n, T neg_inv_n)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T>::is_modulo, "");
    // As in msr_REDC_non_minimized(), protect against undefined behavior:
    using V = typename safely_promote_unsigned<T>::type;
    static_assert(ma::ma_numeric_limits<V>::is_modulo, "");

    static constexpr int bit_width_T = ma::ma_numeric_limits<T>::digits;
    static_assert(bit_width_T % 2 == 0, "");   // bit_width_T divisible by 2
    // MontySqrtRange requires  modulus < sqrt(R)
    static constexpr T sqrtR = static_cast<T>(1) << (bit_width_T / 2);
    HPBC_PRECONDITION2(1 < n && n < sqrtR);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(0 < x && x < sqrtR);
    HPBC_PRECONDITION2(0 < y && y < sqrtR);

    // Since x < sqrtR and y < sqrtR,  x*y < sqrtR*sqrtR == R.
    // We have x*y < R, so x*y will fit in type T without overflow.
    T u_lo = static_cast<T>(static_cast<V>(x) * static_cast<V>(y));
    T result = msr_REDC_non_minimized(u_lo, n, neg_inv_n);

    HPBC_POSTCONDITION2(0 < result && result <= n);
    return result;
}



// MontySqrtRange uses optimizations based on input and output V values being
// 0 < val <= n_, and on modulus < sqrtR.
// These restrictions allow us to implement a more efficient version of the REDC
// algorithm  (in the function msr_REDC_non_minimized()), by omitting some
// conditionals and calculations that would normally be needed.

// The class member variable names are based on the webpage
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
template <typename T>
class MontySqrtRange {
public:
    class MontgomeryValue {
        friend MontySqrtRange;
        explicit MontgomeryValue(T val) : value(val) {}
    public:
#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Weffc++"
#endif
        MontgomeryValue() {} // This constructor purposely does not initialize
        // 'value' - the contents are undefined until the object is assigned to.
#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif
    protected:
        T get() const { return value; }
        T value;
    };
private:
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    const T n_;   // the modulus
    const T r_mod_n_;
    const T neg_inv_n_;
    const T r_squared_mod_n_;
    using V = MontgomeryValue;
public:
    using montvalue_type = V;
    using template_param_type = T;

    explicit MontySqrtRange(T modulus) : n_(modulus),
                             r_mod_n_(getRModN(n_)),
                             neg_inv_n_(negative_inverse_mod_r(n_)),
                             r_squared_mod_n_( modular_arithmetic::
                                    modular_multiplication_prereduced_inputs(
                                                       r_mod_n_, r_mod_n_, n_) )
    {
        namespace ma = hurchalla::modular_arithmetic;
        static constexpr int bitsT = ma::ma_numeric_limits<T>::digits;
        static_assert(bitsT % 2 == 0, "");   // bitsT divisible by 2
        // MontySqrtRange requires  modulus < sqrt(R)
        static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);
        HPBC_PRECONDITION2(1 < modulus && modulus < sqrtR);
        HPBC_PRECONDITION2(modulus % 2 == 1);

        // Note: unityValue == (the montgomery form of 1)==(1*R)%n_ == r_mod_n_.
        // getRModN() guarantees the below.  getUnityValue() and
        // getNegativeOneValue() rely on it.
        HPBC_INVARIANT2(0 < r_mod_n_ && r_mod_n_ < n_);
        // Since n_ == modulus is odd and n_ > 1, n_ can not divide R*R==2^y.
        // Thus  r_squared_mod_n_ == R*R (mod n_) != 0.  convertIn relies on it.
        HPBC_INVARIANT2(0 < r_squared_mod_n_ && r_squared_mod_n_ < n_);
    }
    MontySqrtRange(const MontySqrtRange&) = delete;
    MontySqrtRange& operator=(const MontySqrtRange&) = delete;

    static constexpr T max_modulus()
    {
        static_assert(modular_arithmetic::ma_numeric_limits<T>::digits%2==0,"");
        return (static_cast<T>(1) <<
                      (modular_arithmetic::ma_numeric_limits<T>::digits/2)) - 1;
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
        // Compute R%n.  For example, if R==2^64, arithmetic wraparound behavior
        // of the unsigned integral type T results in (0 - n) representing
        // (2^64 - n).  Thus, rModN = R%n == (2^64)%n == (2^64 - n)%n == (0-n)%n
        T rModN = static_cast<T>(tmp % n);
        // Since n is odd and > 1, n does not divide R==2^x.  Thus, rModN != 0
        HPBC_POSTCONDITION2(0 < rModN && rModN < n);
        return rModN;
    }

public:
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        return (0 < x.get() && x.get() <= n_);
    }

    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        V cfx = getCanonicalValue(x);
        bool good = isValid(x);
        return (x.get() == cfx.get() && good);
    }

    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    // We require a < sqrtR, which is a bit of a hack since MontgomeryForm class
    // expects that any T value >= 0 is ok to use as input for convertIn.
    // Ideally we would address this by renaming this class to something like
    // MontyDoubleWidth, and allowing all MontgomeryValues for this class to be
    // any T value >= 0, while setting
    // T2 = sized_uint<2* ma_numeric_limits<T>::digits>::type, and producing a
    // compile time error if
    // ma_numeric_limits<T2>::digits > HURCHALLA_TARGET_BIT_WIDTH.
    // n_, r_mod_n_, neg_inv_n_, r_squared_mod_n_  would all be type T2, and
    // most function calls in this class would work with type T2 values, and V
    // would wrap type T2.  While this all sounds like a big change, it is not-
    // this class effectively already works like this proposal, with the current
    // type T playing the role of the proposed T2, and convertIn's precondition
    // a < sqrtR playing a psuedo-role of the proposed T.
    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        namespace ma = hurchalla::modular_arithmetic;
        static constexpr int bitsT = ma::ma_numeric_limits<T>::digits;
        static_assert(bitsT % 2 == 0, "");
        static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);
        HPBC_PRECONDITION2(0 <= a && a < sqrtR);

        HPBC_INVARIANT2(1 < n_ && n_ < sqrtR);
        HPBC_INVARIANT2(0 < r_squared_mod_n_ && r_squared_mod_n_ < n_);
        // thus:  0 < r_squared_mod_n_ < sqrtR
        T result;
        HURCHALLA_LIKELY_IF (a > 0) {
            // We have  0 < a < sqrtR  and  0 < r_squared_mod_n_ < sqrtR,  which
            // satisfies the preconditions for msr_montmul_non_minimized().
            result = msr_montmul_non_minimized(a, r_squared_mod_n_, n_,
                                                                    neg_inv_n_);
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

    HURCHALLA_FORCE_INLINE V getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_,
        // and 0 < r_mod_n_ < n_.
        HPBC_POSTCONDITION2(isCanonical(V(r_mod_n_)));
        return V(r_mod_n_);
    }

    HURCHALLA_FORCE_INLINE V getZeroValue() const
    {
        // We want returnVal == (0*R)%n_, but since isValid() requires
        // 0 < returnVal <= n_, we return n_ (n_ ≡ 0 (mod n_))
        HPBC_POSTCONDITION2(isCanonical(V(n_)));
        return V(n_);
    } 

    HURCHALLA_FORCE_INLINE V getNegativeOneValue() const
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

        HPBC_POSTCONDITION2(isCanonical(V(negOne)));
        return V(negOne);
    }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);

        T prod = msr_REDC_non_minimized(x.get(), n_, neg_inv_n_);

        // msr_REDC_non_minimized() postconditions guarantee the following
        HPBC_POSTCONDITION2(0 < prod && prod <= n_);

        T minimized_result;
        HURCHALLA_LIKELY_IF (prod != n_)
            minimized_result = prod;
        else
            minimized_result = static_cast<T>(0);
        HPBC_POSTCONDITION2(minimized_result < n_);
        return minimized_result;
    }

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        return x;
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        HPBC_PRECONDITION2(0 < y.get() && y.get() <= n_);

        // Since we know n < sqrtR  (guaranteed by the constructor),  and x < n
        // and  y < n,  we have  x < sqrtR  and  y < sqrtR,  which satisfies the
        // preconditions of msr_montmul_non_minimized().
        T prod = msr_montmul_non_minimized(x.get(), y.get(), n_, neg_inv_n_);

        // msr_montmul_non_minimized() postconditions guarantee the following
        HPBC_POSTCONDITION2(0 < prod && prod <= n_);
        // Since 0 < prod <= n, we don't want to reduce mod n;  prod is in the
        // canonical form required by most of the class functions.
        return V(prod);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        T a = x.get();
        T b = y.get();
        HPBC_PRECONDITION2(0 < a && a <= n_);
        HPBC_PRECONDITION2(0 < b && b <= n_);
        HPBC_INVARIANT2(n_ > 0);

        T result = montadd_sqrt_range(a, b, n_);

        HPBC_POSTCONDITION2(0 < result && result <= n_);
        return V(result);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        T a = x.get();
        T b = y.get();
        HPBC_PRECONDITION2(0 < a && a <= n_);
        HPBC_PRECONDITION2(0 < b && b <= n_);
        HPBC_INVARIANT2(n_ > 0);

        T result = montsub_sqrt_range(a, b, n_);

        HPBC_POSTCONDITION2(0 < result && result <= n_);
        return V(result);
    }

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        // we can't improve efficiency much over plain subtract,
        // so just delegate to subtract
        return subtract(x, y);
    }

    HURCHALLA_FORCE_INLINE V subtract_canonical_value(V x, V y) const
    {
        // All montgomery values are canonical for this class, so we just
        // delegate to subtract.
        return subtract(x, y);
    }

    HURCHALLA_FORCE_INLINE V add_canonical_value(V x, V y) const
    {
        // All montgomery values are canonical for this class, so we just
        // delegate to add.
        return add(x, y);
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, V z, PTAG) const
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

        V prod = multiply(x, y, PTAG());
        V sum = add(prod, z);

        HPBC_POSTCONDITION2(0 < sum.get() && sum.get() <= n_);
        // Since 0 < sum <= n, we don't want to reduce mod n;  sum is in the
        // canonical form required by most of the class functions.
        return sum;
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, V z, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        HPBC_PRECONDITION2(0 < y.get() && y.get() <= n_);

        // See the optimization discussion inside fmadd() - it applies here too.
        V prod = multiply(x, y, PTAG());
        V diff = subtract(prod, z);

        HPBC_POSTCONDITION2(0 < diff.get() && diff.get() <= n_);
        // Since 0 < diff <= n, we don't want to reduce mod n;  diff is in the
        // canonical form required by most of the class functions.
        return diff;
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V famul(V x, V y, V z, PTAG) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        HPBC_PRECONDITION2(0 < y.get() && y.get() <= n_);
        HPBC_PRECONDITION2(0 < z.get() && z.get() <= n_);

        // Unfortunately for MontySqrtRange, it's not possible to do a simple
        // add of sum=x+y prior to multiplying, since the sum might be greater
        // than sqrtR (e.g. when n is very close to sqrtR and x and n are both
        // very close to n), and that would break a precondition for calling
        // msr_montmul_non_minimized().  Instead we must do a modular addition
        // to get the sum, which means famul() simply wraps this class's add
        // and multiply functions.
        // [Future(?) Note: A hypothetical MontySqrtRangeDiv2 class (requiring
        // modulus < sqrt(R)/2) could provide the optimization of using a simple
        // addition rather than a modular addition: this hypothetical class
        // would have a montmul function that required a*b<R, and using a simple
        // addition  a=(x+y) <= modulus+modulus == 2*modulus, and letting
        // b = z <= modulus, we would have
        // a*b <= 2*modulus*modulus < 2*sqrt(R)*sqrt(R)/4 == R/2, which would
        // satisfy the hypothetical class's montmul requirement.]
        V sum = add(x, y);
        V result = multiply(sum, z, PTAG());

        HPBC_POSTCONDITION2(0 < result.get() && result.get() <= n_);
        return result;
    }
};


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
