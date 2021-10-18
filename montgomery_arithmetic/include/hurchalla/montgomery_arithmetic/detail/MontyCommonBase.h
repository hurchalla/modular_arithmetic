// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/mont_add_canonical_value.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/mont_subtract_canonical_value.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/detail/BaseMontgomeryValue.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// For discussion purposes throughout this file, given an unsigned integral type
// T, let R = 2^(ut_numeric_limits<T>::digits).  For example: if T is uint64_t
// then R = 2^64.  The name 'R' is based on the wikipedia presentation
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
//
// This base class uses the CRTP idiom
// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
// This is the base class shared by most montgomery forms (the experimental
// MontySqrtRange is an exception).
template <template <typename> class Derived, typename T>
class MontyCommonBase {
/*
    MontgomeryValue getMontgomeryValue(WideMontValue x) const
    MontgomeryValue getMontgomeryValue(SquaringValue x) const
    CanonicalValue getCanonicalValue(WideMontValue x) const
    CanonicalValue getCanonicalValue(SquaringValue x) const
    template <class PTAG = LowlatencyTag>
    WideMontValue multiply(WideMontValue x, MontgomeryValue y, bool& resultIsZero) const
    template <class PTAG = LowlatencyTag>
    SquaringValue square(SquaringValue x) const
    template <class PTAG = LowlatencyTag>
    SquaringValue fusedSquareSub(SquaringValue x, CanonicalValue cv) const
    template <class PTAG = LowlatencyTag>
    SquaringValue fusedSquareAdd(SquaringValue x, CanonicalValue cv) const
*/
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    using D = Derived<T>;
 protected:
    struct W : public BaseMontgomeryValue<T> { // wide montgomery value type
        using BaseMontgomeryValue<T>::BaseMontgomeryValue;
    };
    struct SQ : public BaseMontgomeryValue<T> { //squaring montgomery value type
        using BaseMontgomeryValue<T>::BaseMontgomeryValue;
    };
    struct V : public W { using W::W; };  // regular montomery value type
    struct C : public V { using V::V; };  // canonical montgomery value type
    const T n_;   // the modulus
    const T r_mod_n_;
    const T inv_n_;
    const T r_squared_mod_n_;

    explicit MontyCommonBase(T modulus) : n_(modulus),
                       r_mod_n_(get_R_mod_n(n_)),
                       inv_n_(inverse_mod_R(n_)),
                       r_squared_mod_n_(get_Rsquared_mod_n(n_, inv_n_, r_mod_n_,
                                        typename D::MontyTag()))
    {
        HPBC_PRECONDITION2(modulus % 2 == 1);
        HPBC_PRECONDITION2(modulus > 1);
        // Note: unityValue == (the montgomery form of 1)==(1*R)%n_ == r_mod_n_.
        //
        // get_R_mod_n() and get_Rsquared_mod_n() guarantee the below.
        // getUnityValue() and getNegativeOneValue() both rely on it.
        HPBC_INVARIANT2(0 < r_mod_n_ && r_mod_n_ < n_);
        HPBC_INVARIANT2(r_squared_mod_n_ < n_);
    }
    MontyCommonBase(const MontyCommonBase&) = delete;
    MontyCommonBase& operator=(const MontyCommonBase&) = delete;


    // intended for use in preconditions/postconditions
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        T em = static_cast<const D*>(this)->getExtendedModulus();
        return (x.get() < em);
    }

    // intended for use in preconditions/postconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(W x) const
    {
        V cfx = static_cast<const D*>(this)->getCanonicalValue(x);
        // Any fully reduced value (0 <= value < n_) must be canonical.  Class
        // Derived must be implemented to respect this.
        HPBC_INVARIANT2((0 <= x.get() && x.get() < n_) ? x.get() == cfx.get() :
                                                         x.get() != cfx.get());
        return (x.get() == cfx.get());
    }

 public:
    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        HPBC_INVARIANT2(r_squared_mod_n_ < n_);
        // As a precondition, REDC requires  a * r_squared_mod_n < n*R.  This
        // will always be satisfied-  we know from the invariant above that
        // r_squared_mod_n < n.  Since a is a type T variable, we know a < R.
        // Therefore,  a * r_squared_mod_n < n * R.
        T u_lo;
        T u_hi = unsigned_multiply_to_hilo_product(u_lo, a, r_squared_mod_n_);
        // u_hi < n  guarantees we had  a * r_squared_mod_n == u < n*R.  See
        // REDC_non_finalized() in Redc.h for proof.
        HPBC_PRECONDITION2(u_hi < n_);

        T result = REDC(u_hi, u_lo, n_, inv_n_, typename D::MontyTag(),
                                                               LowlatencyTag());
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    HURCHALLA_FORCE_INLINE T convertOut(W x) const
    {
        T u_hi = 0;
        T u_lo = x.get();
        T result= REDC(u_hi, u_lo, n_, inv_n_, FullrangeTag(), LowlatencyTag());
        HPBC_POSTCONDITION2(result < n_);
        return result;
    }

    HURCHALLA_FORCE_INLINE C getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_
        HPBC_INVARIANT2(isCanonical(W(r_mod_n_)));
        return C(r_mod_n_);
    }

    HURCHALLA_FORCE_INLINE C getZeroValue() const
    {
        // zeroValue == (0*R)%n_
        HPBC_INVARIANT2(isCanonical(W(0)));
        return C(0);
    }

    HURCHALLA_FORCE_INLINE C getNegativeOneValue() const
    {
        // We want to get returnVal = getCanonicalValue(subtract(getZeroValue(),
        //                                               getUnityValue())).
        //   getZeroValue() returns a value belonging to the equivalence class
        //   0*R (mod n_).  This equivalence class can equally be represented by
        //   the value n_ (mod n_).  getUnityValue() returns a value belonging
        //   to the equivalence class 1*R (mod n_), which is the same as 
        //   r_mod_n_ (mod n_).  Therefore the subtraction results in the
        //   equivalence class (n_ - r_mod_n_) (mod n_).
        //   The constructor established the invariant  0 < r_mod_n_ < n_
        //   Thus we also know  0 < n_ - r_mod_n_ < n_.  This means
        //   (n_ - r_mod_n_)  is fully reduced, and thus canonical.
        HPBC_INVARIANT2(0 < r_mod_n_ && r_mod_n_ < n_);
        T ret = static_cast<T>(n_ - r_mod_n_);
        HPBC_ASSERT2(0 < ret && ret < n_);

        HPBC_POSTCONDITION2(isCanonical(W(ret)));
        return C(ret);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T em = static_cast<const D*>(this)->getExtendedModulus();
        T z = modular_addition_prereduced_inputs(x.get(), y.get(), em);
        HPBC_POSTCONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V add_canonical_value(V x, C y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isCanonical(y));
        HPBC_PRECONDITION2(y.get() < n_); // isCanonical() should guarantee this
        T z = mont_add_canonical_value<T>::call(x.get(), y.get(), n_);
        // mont_add_canonical_value guarantees that z <= std::max(x, n-1).
        // Thus if x < n, then z < n.  Or in other words, if x is canonical,
        // then z is canonical.
        HPBC_POSTCONDITION2(isCanonical(x) ? isCanonical(V(z)) : true);
        HPBC_POSTCONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T em = static_cast<const D*>(this)->getExtendedModulus();
        T z = modular_subtraction_prereduced_inputs(x.get(), y.get(), em);
        HPBC_POSTCONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V subtract_canonical_value(V x, C y) const
    {
        HPBC_PRECONDITION2(isCanonical(y));
        HPBC_PRECONDITION2(y.get() < n_); // isCanonical() should guarantee this
        HPBC_PRECONDITION2(isValid(x));
        T z = mont_subtract_canonical_value<T>::call(x.get(), y.get(), n_);
        // mont_subtract_canonical_value guarantees that z <= std::max(x, n-1).
        // Thus if x < n, then z < n.  Or in other words, if x is canonical,
        // then z is canonical.
        HPBC_POSTCONDITION2(isCanonical(x) ? isCanonical(V(z)) : true);
        HPBC_POSTCONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE C subtract_dual_canonical_values(C x, C y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(x.get() < n_); // isCanonical() should guarantee this
        HPBC_PRECONDITION2(isCanonical(y));
        HPBC_PRECONDITION2(y.get() < n_); // isCanonical() should guarantee this
        T z = modular_subtraction_prereduced_inputs(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(isCanonical(W(z)));
        return C(z);
    }

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T result = absolute_value_difference(x.get(), y.get());
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    // Multiplies two montgomery values x and y, and returns the product as a
    // montgomery value.  Sets isZero to indicate if the product ≡ 0 (mod n).
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, bool& isZero, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        // As a precondition, REDC requires  x*y < n*R.  This will always be
        // satisfied for all Monty classes known to derive from this class --
        // For MontyFullRange: its constructor requires modulus < R, which
        //   means n < R.  For MontyFullRange, isValid(a) returns (a < n), so
        //   we know by this function's preconditions that x < n and y < n.
        //   Therefore  x*y < n*n < n*R.
        // For MontyHalfRange: its constructor requires modulus < R/2, which
        //   means n < R/2 < R.  For MontyHalfRange, isValid(a) returns (a < n),
        //   so we know by this function's preconditions that x < n and y < n.
        //   Therefore  x*y < n*n < n*R.
        // For MontyQuarterRange: its constructor requires modulus < R/4, which
        //   means n < R/4.  For MontyQuarterRange isValid(a) returns (a < 2*n),
        //   so we know by this function's preconditions that x < 2*n and
        //   y < 2*n.  Thus  x*y < (2*n)*(2*n) == 4*n*n < 4*n*R/4 == n*R.
        T u_lo;
        T u_hi = unsigned_multiply_to_hilo_product(u_lo, x.get(), y.get());
        // u_hi < n  implies that  x*y == u < n*R.  See REDC_non_finalized()
        // in Redc.h for proof.
        HPBC_ASSERT2(u_hi < n_);

        T result = REDC(u_hi, u_lo, n_, inv_n_, isZero, typename D::MontyTag(),
                                                                        PTAG());
        HPBC_POSTCONDITION2(isZero ==
             (static_cast<const D*>(this)->getCanonicalValue(V(result)).get() ==
              getZeroValue().get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    // Multiplies two montgomery values x and y, and then subtracts canonical
    // value z from the product.  Returns the resulting montgomery value.
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, C z, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_PRECONDITION2(isCanonical(z));
        HPBC_PRECONDITION2(z.get() < n_); // isCanonical() should guarantee this
        T u_lo;
        T u_hi = unsigned_multiply_to_hilo_product(u_lo, x.get(), y.get());
        // Assuming theoretical unlimited precision standard multiplication,
        // REDC requires  u = x*y < n*R.  See multiply() for why this function
        // will always satisfy the requirement.  u_hi < n guarantees we had
        // x*y == u < n*R.  See REDC_non_finalized() in Redc.h for proof.
        HPBC_ASSERT2(u_hi < n_);

        // Performing the modular sub prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // subtraction.  See fmadd() below for a proof which could be easily
        // adapted for the modular subtraction in fmsub().
        // The following calculations should execute in parallel with the first
        // two multiplies in REDC(), since those mutiplies do not depend on
        // these calculations.  (Instruction level parallelism)
        T diff = modular_subtraction_prereduced_inputs(u_hi, z.get(), n_);
        HPBC_ASSERT2(diff < n_);
        T result = REDC(diff, u_lo, n_, inv_n_, typename D::MontyTag(), PTAG());

        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    // Multiplies two montgomery values x and y, and then adds canonical
    // value z to the product.  Returns the resulting montgomery value.
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, C z, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_PRECONDITION2(isCanonical(z));
        HPBC_PRECONDITION2(z.get() < n_); // isCanonical() should guarantee this
        T u_lo;
        T u_hi = unsigned_multiply_to_hilo_product(u_lo, x.get(), y.get());
        // Assuming theoretical unlimited precision standard multiplication,
        // REDC requires  u = x*y < n*R.  See multiply() for why this function
        // will always satisfy the requirement.  u_hi < n guarantees we had
        // x*y == u < n*R.  See REDC_non_finalized() in Redc.h for proof.
        HPBC_ASSERT2(u_hi < n_);

        // The most obvious way to carry out this function would be as follows:
        // result = REDC(u_hi, u_lo, n_, inv_n_, typename D::MontyTag(), PTAG())
        // and then performing a modular addition of result with z
        // sum = this->add(result, z).
        // However we'll show we can perform a modular addition first, and then
        // the REDC, and our final result will be congruent to sum above, and
        // our final result will satisfy
        // this->isValid(final_result) == true.  If we can achieve this, then
        // for our purposes final_result is equivalent to the sum above.
        //
        // The advantage of this alternate method is that the modular addition
        // should execute in parallel with the first two multiplies inside
        // REDC(), because the mutiplies do not depend on the result of the
        // addition (i.e. instruction level parallelism).  In contrast, the most
        // obvious way of implementing this function, which we described above,
        // can not take advantage of instruction level parallelism because the
        // inputs to the ending modular addition depend upon the result of the
        // REDC.
        //
        // Let's prove final_result is congruent to sum and that it is valid.
        // Proof:
        // Let Rinverse ≡ R^(-1) (mod n_)       ("R^(-1)" is R to the power -1)
        // Note that since R is a power of 2 and n_ is odd, Rinverse always
        // exists.  Let u = u_hi*R + u_lo.  If we call
        // result = REDC(u_hi, u_lo, n_, inv_n_, typename D::MontyTag(), PTAG())
        // REDC guarantees that result satisfies
        // result     ≡ u*Rinverse  (mod n_).   And therefore,
        // result     ≡ (u_hi*R + u_lo)*Rinverse  (mod n_)
        // result + z ≡ (u_hi*R + u_lo)*Rinverse + z  (mod n_)
        // result + z ≡ (u_hi*R + u_lo)*Rinverse + z*R*Rinverse  (mod n_)
        // result + z ≡ ((u_hi + z)*R + u_lo)*Rinverse  (mod n_)
        // result + z ≡ (((u_hi + z)%n_)*R + u_lo)*Rinverse  (mod n_)
        // The earlier  sum = this->add(result, z)  satisfies
        // sum ≡ result + z  (mod n_),  and so
        // sum        ≡ (((u_hi + z)%n_)*R + u_lo)*Rinverse  (mod n_)
        //
        // Earlier in this function, we showed by precondition and assertion
        // that  u_hi < n_, and z < n_.  Therefore we can perform (u_hi + z)%n_
        // via the function modular_addition_prereduced_inputs(), as follows:
        T v_hi = modular_addition_prereduced_inputs(u_hi, z.get(), n_);
        // By substitution we have 
        // sum        ≡ (v_hi*R + u_lo)*Rinverse  (mod n_)
        //
        // Because v_hi == (u_hi + z)%n_,  we know  v_hi < n_.
        HPBC_ASSERT2(v_hi < n_);
        // Thus  v_hi <= n_ - 1,  and  v_hi*R <= n_*R - R.  And
        // v_hi*R + u_lo <= n_*R - R + u_lo.  Since u_lo is type T, it
        // satisfies u_lo < R.  And thus
        // v_hi*R + u_lo <= n_*R - R + u_lo < n_*R - R + R == n_*R.
        // Therefore we know  (v_hi*R + u_lo) < n_*R.
        //
        // Let  v == v_hi*R + u_lo.  We saw just above that v < n_*R, which
        // satisfies REDC's precondition requiring an input < n_*R.  Therefore
        // we can call REDC with the input v == (v_hi*R + u_lo), and the REDC
        // algorithm guarantees it will return the value
        // final_result ≡ (v_hi*R + u_lo)*Rinverse  (mod n_).  And thus,
        // final_result ≡ sum  (mod n_).
        T final_result =
                   REDC(v_hi, u_lo, n_, inv_n_, typename D::MontyTag(), PTAG());
        // REDC's postcondition guarantees it will return a valid value for
        // the given MontyTag, and so
        HPBC_POSTCONDITION2(isValid(V(final_result)));
        // We showed  final_result ≡ sum  (mod n_),  and that
        // final_result is valid.  For our purposes this means final_result is
        // equivalent to sum, and thus we can return it instead of sum.
        return V(final_result);
    }

    // Returns the greatest common divisor of the standard representations
    // (non-montgomery) of both x and the modulus, using the supplied functor.
    // The functor must take two integral arguments of the same type and return
    // the gcd of those two arguments.
    template <class F>
    HURCHALLA_FORCE_INLINE T gcd_with_modulus(W x, const F& gcd_functor) const
    {
        // Proof that gcd(x.get(), n_) == gcd(convertOut(x), n_)
        // -----------------------------------------------------
        // Let the integer g = x.get(), and let the integer c = convertOut(x).
        // Let the integer d be a divisor of n_.
        // We will use mathematical integers (with infinite precision and no
        // overflow) throughout this discussion.
        //
        // Because g is a value in montgomery domain, we know g ≡ c*R (mod n_),
        // and thus there exists some integer k such that  g == c*R + k*n_.
        // Since n_ (by constructor precondition) is odd, we know n_ and R are
        // coprime, and thus d can not be a divisor of R (unless d==1, in which
        // case d divides all integers).
        // Therefore d divides c*R if and only if d divides c.
        // Assume d divides c:
        //    Then d divides c*R.  Since d divides n_, d also divides k*n_.
        //    Thus d divides c*R + k*n_ == g.
        // Assume d divides g:
        //    Since d divides n_, d also divides k*n_.  Thus d divides
        //    g - k*n_ == c*R.  Since we showed d can not be a divisor of R
        //    (unless d==1), d must divide c.
        // Therefore d divides g if and only if d divides c.
        //
        // Let p be the greatest common divisor of g and n_.  Since p
        // divides g, p must divide c.  Let  q = gcd(c, n_).  Since q
        // divides c, q must divide g.  Since q also divides n_, q is a
        // common divisor of g and n_, and thus q <= p.  Since p divides
        // both c and n_, p is a common divisor of c and n_, and thus p <= q.
        // Since q <= p and p <= q, we know q == p.  Therefore:
        // gcd(g, n_) == gcd(c, n_).
        // -----------------------------------------------
        //
        // We want to return the value  q = gcd(convertOut(x), n_).
        // By the proof above, we can instead return the equivalent value
        // p = gcd(x.get(), n_)  which we can compute more efficiently.
        T p = gcd_functor(x.get(), n_);
        // Our postconditions assume the functor implementation is correct.
        HPBC_POSTCONDITION2(0 < p && p <= n_ && (x.get() == 0 || p <= x.get()));
        HPBC_POSTCONDITION2(n_ % p == 0);
        HPBC_POSTCONDITION2(x.get() % p == 0);
        return p;
    }
};


}} // end namespace

#endif
