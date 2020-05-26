
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SQRT_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_SQRT_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/montgomery_arithmetic/detail/make_safe_unsigned_integer.h"
#include "hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h"
#include "hurchalla/montgomery_arithmetic/detail/MontgomeryValue.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// This function is based on REDC_non_minimized() from monty_common.h.  It is
// altered to omit calculations that are not needed, given the preconditions of
// n < sqrt(R), and 0 < x <= n, and 0 < y <= n.
template <typename T>
HURCHALLA_FORCE_INLINE T msr_montmul_non_minimized(T x, T y, T n, T neg_inv_n)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");

    // For casts, we want to use types that are protected from surprises and
    // undefined behavior due to the unsigned integral promotion rules in C++.
    // https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/
    using V = typename make_safe_unsigned_integer<T>::type;
    static_assert(std::numeric_limits<V>::is_modulo, "");

    static constexpr int bit_width_T = std::numeric_limits<T>::digits;
    static_assert(bit_width_T % 2 == 0, "");   // bit_width_T divisible by 2
    // MontySqrtRange requires  modulus < sqrt(R)
    static constexpr T sqrtR = static_cast<T>(1) << (bit_width_T / 2);
    HPBC_PRECONDITION2(1 < n && n < sqrtR);
    HPBC_PRECONDITION2(n % 2 == 1);
    HPBC_PRECONDITION2(0 < x && x <= n);
    HPBC_PRECONDITION2(0 < y && y <= n);

    // assert(n * neg_inv_n ≡ -1 (mod R))
    HPBC_PRECONDITION2(
                static_cast<T>(static_cast<V>(n) * static_cast<V>(neg_inv_n)) ==
                static_cast<T>(static_cast<V>(0) - static_cast<V>(1))
                );

    // Since n < sqrtR, and x <= n and y <= n,  x*y <= n*n < sqrtR*sqrtR == R.
    // We have x*y < R, and since n>1, x*y < R < R*n.  Thus we've satisfied the
    // basic requirement for montgomery multiplication that u = x*y < n*R.
    // Since u = x*y < R, we have  u_lo = (x*y)%R == x*y, and u_hi == 0.
    V u_lo = static_cast<V>(x) * static_cast<V>(y);

    // compute  m = (u * neg_inv_n) % R
    T m = static_cast<T>(u_lo * static_cast<V>(neg_inv_n));

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
    // However, we know  u_lo = (x*y)%R, and we proved  u_lo == x*y < R.  Since
    // our preconditions specify x>0 and y>0, we know  x*y > 0.  Thus  u_lo > 0,
    // or more specifically, u_lo != 0.  The calculation of t_hi simplifies to
    t_hi += static_cast<T>(1);

    // REDC_non_minimized() would normally next calculate
    // ovf = (t_hi < u_hi);
    // But we know u_hi == 0, so ovf = (t_hi < u_hi) == (t_hi < 0) == false.
    // Thus  ovf = false.

    // The discussion prior to Assertion #1 proves that mn_hi < n, and therefore
    // 0 < mn_hi + 1 < n + 1.  Since t_hi = mn_hi + 1, we know  0 < t_hi <= n.
    HPBC_POSTCONDITION2(0 < t_hi && t_hi <= n);
    // From REDC_non_minimized() we have the postcondition:
    //   T minimized_result = (ovf || t_hi >= n) ? (t_hi - n) : t_hi;
    //   HPBC_POSTCONDITION2(minimized_result < n);
    // since  ovf == false  and  0 < t_hi <= n,  we can simplify this to
    if (HPBC_POSTCONDITION2_MACRO_IS_ACTIVE) {
        T minimized_result = (t_hi == n) ? 0 : t_hi;
        HPBC_POSTCONDITION2(minimized_result < n);
    }

    // return the non-minimized result
    return t_hi;
}



// MontySqrtRange uses optimizations based on input and output V values being
// 0 < val <= n_, and on modulus < sqrtR.
// These restrictions allow us to implement a more efficient version of the REDC
// algorithm  (in the function msr_montmul_non_minimized()), by omitting some
// conditionals and calculations that would normally be needed.

// The class member variable names are based on the webpage
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
//
// For discussion purposes, let R = 2^(std::numeric_limits<T>::digits).  For
// example if T is uint64_t, then R = 2^64.
template <typename T>
class MontySqrtRange final {
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    const T n_;   // the modulus
    const T neg_inv_n_;
    const T r_mod_n_;
    const T r_squared_mod_n_;
public:
    using V = MontgomeryValue<T>;
    using montvalue_type = V;
    using template_param_type = T;

    explicit MontySqrtRange(T modulus) : n_(modulus),
                             neg_inv_n_(negative_inverse_mod_r(n_)),
                             r_mod_n_(getRModN(n_)),
                             r_squared_mod_n_( modular_arithmetic::
                                    modular_multiplication_prereduced_inputs(
                                                       r_mod_n_, r_mod_n_, n_) )
    {
        static constexpr int bitsT = std::numeric_limits<T>::digits;
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
        constexpr int bitsT = std::numeric_limits<T>::digits;
        static_assert(bitsT % 2 == 0, "");   // bitsT divisible by 2
        constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);
        return sqrtR - 1;
    }

private:
    static T getRModN(T n)
    {
        HPBC_PRECONDITION2(n % 2 == 1);
        HPBC_PRECONDITION2(n > 1);
        // Assign a tmp T variable rather than directly using the intermediate
        // expression, in order to avoid a negative value (and a wrong answer)
        // in cases where 'n' would be promoted to type 'int'.
        T tmp = static_cast<T>(0) - n;
        // Compute R%n.  For example, if R==2^64, arithmetic wraparound behavior
        // of the unsigned integral type T results in (0 - n) representing
        // (2^64 - n).  Thus, rModN = R%n == (2^64)%n == (2^64 - n)%n == (0-n)%n
        T rModN = tmp % n;
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
        V cfx = getCanonicalForm(x);
        bool good = isValid(x);
        return (x == cfx && good); 
    }

    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        HPBC_INVARIANT2(0 < r_squared_mod_n_ && r_squared_mod_n_ < n_);
        HPBC_PRECONDITION2(0 <= a && a < n_);
        // multiply requires valid input values, and 0 is the single possible
        // invalid value of 'a' for the multiply.  We treat this case a == 0
        // separately, with  a*R (mod n) ≡ 0*R (mod n) ≡ 0 (mod n) ≡ n (mod n).
        V result = (a > 0) ? multiply(V(a), V(r_squared_mod_n_)) : V(n_);
        HPBC_POSTCONDITION2(0 < result.get() && result.get() <= n_);
        return result;
    }

    HURCHALLA_FORCE_INLINE V getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_,
        // and 0 < r_mod_n_ < n_.
        HPBC_INVARIANT2(isCanonical(V(r_mod_n_)));
        return V(r_mod_n_);
    }

    HURCHALLA_FORCE_INLINE V getZeroValue() const
    {
        // We want returnVal == (0*R)%n_, but since isValid() requires
        // 0 < returnVal <= n_, we return n_ (n_ ≡ 0 (mod n_))
        V zero(n_);
        HPBC_INVARIANT2(isCanonical(zero));
        return zero;
    } 

    HURCHALLA_FORCE_INLINE V getNegativeOneValue() const
    {
        // We want to get  returnVal = getCanonicalForm(subtract(getZeroValue(),
        //                                               getUnityValue())).
        //   getZeroValue() returns n_, and getUnityValue() returns r_mod_n_.
        //   Therefore the subtraction results in the equivalence class
        //   (n_ - r_mod_n_) (mod n_). The constructor established the invariant
        //   0 < r_mod_n_ < n_.  Thus we know  0 < n_ - r_mod_n_ < n_.  This
        //   means (n_ - r_mod_n_)  satisfies isValid() and getCanonicalForm().
        HPBC_INVARIANT2(n_ > r_mod_n_);
        T negOne = n_ - r_mod_n_;
        HPBC_ASSERT2(0 < negOne && negOne < n_);
        HPBC_INVARIANT2(isCanonical(V(negOne)));
        return V(negOne);
    }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);

        T y = static_cast<T>(1);
        T prod = msr_montmul_non_minimized(x.get(), y, n_, neg_inv_n_);

        // msr_montmul_non_minimized() postconditions guarantee the following
        HPBC_POSTCONDITION2(0 < prod && prod <= n_);
        T minimized_result = (prod != n_) ? prod : 0;
        HPBC_POSTCONDITION2(minimized_result < n_);
        return minimized_result;
    }

    HURCHALLA_FORCE_INLINE V getCanonicalForm(V x) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        return x;
    }

    HURCHALLA_FORCE_INLINE V multiply(V x, V y) const
    {
        HPBC_PRECONDITION2(0 < x.get() && x.get() <= n_);
        HPBC_PRECONDITION2(0 < y.get() && y.get() <= n_);

        T prod = msr_montmul_non_minimized(x.get(), y.get(), n_, neg_inv_n_);

        // msr_montmul_non_minimized() postconditions guarantee the following
        HPBC_POSTCONDITION2(0 < prod && prod <= n_);
        // Since 0 < prod <= n, we don't want to reduce mod n;  prod is in the
        // canonical form required by most of the class functions.
        return V(prod);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        // modular addition (mod n_), except that a result of 0 becomes n_.
        // This is adapted from  modular_addition_prereduced_inputs():
        T a = x.get();
        T b = y.get();
        HPBC_PRECONDITION2(0 < a && a <= n_);
        HPBC_PRECONDITION2(0 < b && b <= n_);
        HPBC_INVARIANT2(n_ > 0);

        T tmp = n_ - b;
        T result = (a <= tmp) ? a+b : a-tmp;

        HPBC_POSTCONDITION2(0 < result && result <= n_);
        return V(result);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        // modular subtraction (mod n_), except that a result of 0 becomes n_.
        // This is adapted from  modular_subtraction_prereduced_inputs():
        T a = x.get();
        T b = y.get();
        HPBC_PRECONDITION2(0 < a && a <= n_);
        HPBC_PRECONDITION2(0 < b && b <= n_);
        HPBC_INVARIANT2(n_ > 0);

        T result = (a>b) ? a-b : n_ - (b-a);

        HPBC_POSTCONDITION2(0 < result && result <= n_);
        return V(result);
    }
};


}} // end namespace

#endif
