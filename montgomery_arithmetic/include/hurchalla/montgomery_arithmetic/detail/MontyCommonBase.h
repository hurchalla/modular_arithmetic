// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/detail/BaseMontgomeryValue.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// For discussion purposes throughout this file, given an unsigned integral type
// T, let R = 1<<(ut_numeric_limits<T>::digits).  For example: if T is uint64_t
// then R = 1<<64.  The name 'R' is based on the wikipedia presentation
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
//
// This base class uses the CRTP idiom
// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
// This is the base class shared by most montgomery forms (the experimental
// MontySqrtRange is an exception).
template <template<typename> class Derived,
          template<typename> class MVTypes,
          typename T>
class MontyCommonBase {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    using D = Derived<T>;
 protected:
    using V = typename MVTypes<T>::V;   // the MontgomeryValue type
    using C = typename MVTypes<T>::C;   // the CanonicalValue type
    const T n_;   // the modulus
    const T r_mod_n_;
    const T inv_n_;
    const T r_squared_mod_n_;

    explicit MontyCommonBase(T modulus) :
         n_(modulus),
         r_mod_n_(::hurchalla::get_R_mod_n(n_)),
         inv_n_(::hurchalla::inverse_mod_R(n_)),
         r_squared_mod_n_(::hurchalla::get_Rsquared_mod_n(n_,inv_n_,r_mod_n_))
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

 public:
    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        HPBC_INVARIANT2(r_squared_mod_n_ < n_);
        // As a precondition, REDC requires  a * r_squared_mod_n < n*R.  This
        // will always be satisfied-  we know from the invariant above that
        // r_squared_mod_n < n.  Since a is a type T variable, we know a < R.
        // Therefore,  a * r_squared_mod_n < n * a < n * R.
        T u_lo;
        T u_hi = ::hurchalla::unsigned_multiply_to_hilo_product(
                                                     u_lo, a, r_squared_mod_n_);
        // Let u = a * r_squared_mod_n.  When u_hi < n, we always have u < n*R.
        // See RedcIncomplete() in ImplRedc.h for proof.
        HPBC_ASSERT2(u_hi < n_);
        const D* child = static_cast<const D*>(this);
        V result = child->montyREDC(u_hi, u_lo, LowlatencyTag());
        HPBC_POSTCONDITION2(child->isValid(result));
        return result;
    }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        T u_hi = 0;
        // get a Natural number (i.e. number >= 0) congruent to x (mod n)
        T u_lo = static_cast<const D*>(this)->getNaturalEquivalence(x);
        namespace hc = ::hurchalla;
        T result = hc::REDC_standard(u_hi, u_lo, n_, inv_n_, LowlatencyTag());

        HPBC_POSTCONDITION2(result < n_);
        return result;
    }

    HURCHALLA_FORCE_INLINE C getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_
        HPBC_INVARIANT2(r_mod_n_ < n_);
        return C(r_mod_n_);
    }

    HURCHALLA_FORCE_INLINE C getZeroValue() const
    {
        return C(0);  // zeroValue == (0*R)%n_
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
        return C(ret);
    }


    template <class PTAG>
    HURCHALLA_FORCE_INLINE V square(V x, PTAG) const
    {
        const D* child = static_cast<const D*>(this);
        HPBC_PRECONDITION2(child->isValid(x));
        T u_lo;
        T u_hi = child->squareToHiLo(u_lo, x);
        HPBC_ASSERT2(u_hi < n_);  // verifies  u < n*R, as required for REDC
        V result = child->montyREDC(u_hi, u_lo, PTAG());
        HPBC_POSTCONDITION2(child->isValid(result));
        return result;
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareSub(V x, C cv, PTAG) const
    {
        const D* child = static_cast<const D*>(this);
        HPBC_PRECONDITION2(child->isValid(x));
        HPBC_INVARIANT2(cv.get() < n_);
        T u_lo;
        T u_hi = child->squareToHiLo(u_lo, x);
        HPBC_ASSERT2(0 <= u_hi && u_hi < n_);
        // Performing the modular sub prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // subtraction.  See fmadd() below for a proof which could be easily
        // adapted for the modular subtraction in here; it also describes
        // why this order of operations is more efficient.
#if 0
        namespace hc = ::hurchalla;
        u_hi = hc::modular_subtraction_prereduced_inputs(u_hi, cv.get(), n_);
#else
        // this provides the same results as the #if code, but can be more
        // efficient for some monty types
        HPBC_ASSERT2(child->isCanonical(C(u_hi)));
        u_hi = child->subtract(C(u_hi), cv).get();
#endif
        HPBC_ASSERT2(u_hi < n_);  // verifies  u < n*R, as required for REDC
        V result = child->montyREDC(u_hi, u_lo, PTAG());
        HPBC_POSTCONDITION2(child->isValid(result));
        return result;
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareAdd(V x, C cv, PTAG) const
    {
        const D* child = static_cast<const D*>(this);
        HPBC_PRECONDITION2(child->isValid(x));
        HPBC_INVARIANT2(cv.get() < n_);
        T u_lo;
        T u_hi = child->squareToHiLo(u_lo, x);
        HPBC_ASSERT2(0 <= u_hi && u_hi < n_);
        // Performing the modular add prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // addition.  See fmadd() below for a proof which applies directly
        // to this function; it also describes why this order of operations
        // is more efficient.
#if 0
        namespace hc = ::hurchalla;
        u_hi = hc::modular_addition_prereduced_inputs(u_hi, cv.get(), n_);
#else
        // this provides the same results as the #if code, but can be more
        // efficient for some monty types
        HPBC_ASSERT2(child->isCanonical(C(u_hi)));
        u_hi = child->add(C(u_hi), cv).get();
#endif
        HPBC_ASSERT2(u_hi < n_);  // verifies  u < n*R, as required for REDC
        V result = child->montyREDC(u_hi, u_lo, PTAG());
        HPBC_POSTCONDITION2(child->isValid(result));
        return result;
    }


    // Multiplies two montgomery values x and y, and returns the product as a
    // montgomery value.  Sets isZero to indicate if the product ≡ 0 (mod n).
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, bool& isZero, PTAG) const
    {
        const D* child = static_cast<const D*>(this);
        HPBC_PRECONDITION2(child->isValid(x));
        HPBC_PRECONDITION2(child->isValid(y));
        T u_lo;
        T u_hi = child->multiplyToHiLo(u_lo, x, y);
        // u_hi < n  implies that  x*y == u < n*R.  See RedcIncomplete() in
        // ImplRedc.h for proof.
        HPBC_ASSERT2(u_hi < n_);  // verify u < n*R, as required for REDC
        V result = child->montyREDC(isZero, u_hi, u_lo, PTAG());
        HPBC_POSTCONDITION2(child->isValid(result));
        HPBC_POSTCONDITION2(isZero ==
                          (child->getCanonicalValue(result) == getZeroValue()));
        return result;
    }

// the Derived class will have overloads of fmsub and fmadd that might offer
// better efficiency, using a dedicated type FV (FusingValue) rather than C
// (CanonicalValue) for the addend/subtrahend z.

    // Multiplies two montgomery values x and y, and then subtracts canonical
    // value z from the product.  Returns the resulting montgomery value.
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, C z, PTAG) const
    {
        const D* child = static_cast<const D*>(this);
        HPBC_PRECONDITION2(child->isValid(x));
        HPBC_PRECONDITION2(child->isValid(y));
        HPBC_INVARIANT2(z.get() < n_);
        T u_lo;
        T u_hi = child->multiplyToHiLo(u_lo, x, y);
        // Assuming theoretical unlimited precision standard multiplication,
        // REDC requires  u = x*y < n*R.  multiplyToHiLo() guarantees that we
        // will receive values for  u_hi, u_lo  that satisfy this requirement.
        // u_hi < n implies that  x*y == u < n*R.  See RedcIncomplete() in
        // ImplRedc.h for proof.
        HPBC_ASSERT2(0 <= u_hi && u_hi < n_);

        // Performing the modular sub prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // subtraction.  See fmadd() below for a proof which could be easily
        // adapted for the modular subtraction in here.
        // Note that the modular subtraction will usually execute in parallel
        // with the first two multiplies in REDC(), since those mutiplies do not
        // depend on the subtraction result.  (Instruction level parallelism)
#if 0
        namespace hc = ::hurchalla;
        u_hi = hc::modular_subtraction_prereduced_inputs(u_hi, z.get(), n_);
#else
        // this provides the same results as the #if code, but can be more
        // efficient for some monty types
        HPBC_ASSERT2(child->isCanonical(C(u_hi)));
        u_hi = child->subtract(C(u_hi), z).get();
#endif
        HPBC_ASSERT2(u_hi < n_);
        V result = child->montyREDC(u_hi, u_lo, PTAG());

        HPBC_POSTCONDITION2(child->isValid(result));
        return result;
    }
    // Multiplies two montgomery values x and y, and then adds canonical
    // value z to the product.  Returns the resulting montgomery value.
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, C z, PTAG) const
    {
        const D* child = static_cast<const D*>(this);
        HPBC_PRECONDITION2(child->isValid(x));
        HPBC_PRECONDITION2(child->isValid(y));
        HPBC_INVARIANT2(z.get() < n_);
        T u_lo;
        T u_hi = child->multiplyToHiLo(u_lo, x, y);
        // Assuming theoretical unlimited precision standard multiplication,
        // REDC requires  u = u_hi*R + u_lo = x*y < n*R.  multiplyToHiLo()
        // guarantees that  u_hi and u_lo  will satisfy this requirement.
        // u_hi < n implies that  x*y == u < n*R.  See RedcIncomplete() in
        // ImplRedc.h for proof.
        HPBC_ASSERT2(0 <= u_hi && u_hi < n_);

        // The most obvious way to carry out this function would be as follows:
        // result = this->montyREDC(u_hi, u_lo, PTAG())
        // and then performing a modular addition of result with z
        // sum = this->add(result, z).
        // However we'll show we can perform a modular addition first, and then
        // the REDC, and our final result will be congruent to sum above, and
        // our final result will satisfy
        // isValid(final_result) == true.  If we can achieve this, then for our
        // purposes final_result is equivalent to the sum above.
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
        // exists.  Let u = u_hi*R + u_lo.  If we let
        // result = this->montyREDC(u_hi, u_lo, PTAG())
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
#if 0
        namespace hc = ::hurchalla;
        T v_hi = hc::modular_addition_prereduced_inputs(u_hi, z.get(), n_);
#else
        // this provides the same results as the #if code, but can be more
        // efficient for some monty types
        HPBC_ASSERT2(child->isCanonical(C(u_hi)));
        T v_hi = child->add(C(u_hi), z).get();
#endif
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
        V final_result = child->montyREDC(v_hi, u_lo, PTAG());
        // REDC's postcondition guarantees it will return a valid value, and so
        HPBC_POSTCONDITION2(child->isValid(final_result));
        // We showed  final_result ≡ sum  (mod n_),  and that
        // final_result is valid.  For our purposes this means final_result is
        // equivalent to sum, and thus we can return it instead of sum.
        return final_result;
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
        // Proof that GCD(getNaturalEquivalence(x), n) == GCD(convertOut(x), n)
        // --------------------------------------------------------------------
        // Let the integer g = ((const D*)this)->getNaturalEquivalence(x),
        // and let the integer c = convertOut(x).
        // Let the integer d be a divisor of n_.
        // We will use mathematical integers (with infinite precision and no
        // overflow) throughout this discussion.
        //
        // GetNaturalEquivalence(x) guarantees that its returned integer g
        // satisfies  g ≡ x (mod n_), and because x is a value in Montgomery
        // domain, we know g ≡ x ≡ c*R (mod n_).  Thus there exists some integer
        // k such that  g == c*R + k*n_.
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
        // p = gcd(((const D*)this)->getNaturalEquivalence(x), n_)
        // which we can compute more efficiently.
        T g = static_cast<const D*>(this)->getNaturalEquivalence(x);
        T p = gcd_functor(g, n_);
        // Our postconditions assume the functor implementation is correct.
        HPBC_POSTCONDITION2(0 < p && p <= n_ && (g == 0 || p <= g));
        HPBC_POSTCONDITION2(n_ % p == 0);
        HPBC_POSTCONDITION2(g % p == 0);
        return p;
    }
};


}} // end namespace

#endif
