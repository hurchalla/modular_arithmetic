// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_HALF_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_HALF_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/halfrange_get_canonical.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/signed_multiply_to_hilo_product.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>


namespace hurchalla { namespace detail {


// The name "Halfrange" signifies that the modulus must be less than R/2, where
// R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t then
// R = 1<<64 and R/2 == 1<<63, and thus it would require  modulus < (1<<63).


struct TagMontyHalfrange final {};


// struct used internally by MontyHalfRange
template <typename T>
struct MontyHRValueTypes {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    using SignedT = typename extensible_make_signed<T>::type;
    static_assert(static_cast<SignedT>(static_cast<T>(static_cast<SignedT>(-1)))
                  == static_cast<SignedT>(-1), "Casting a signed value to "
                   "unsigned and back again must result in the original value");
    static_assert(static_cast<SignedT>(-1) == ~(static_cast<SignedT>(0)),
                            "SignedT must use two's complement representation");
    struct C; // forward declare C so that V can friend it
    // regular montgomery value type
    struct V : public BaseMontgomeryValue<SignedT> {
        HURCHALLA_FORCE_INLINE V() = default;
     protected:
        friend struct C;
        template <typename> friend class MontyHalfRange;
        HURCHALLA_FORCE_INLINE
        explicit V(SignedT a) : BaseMontgomeryValue<SignedT>(a) {}
    };
    // canonical montgomery value type
    struct C : public BaseMontgomeryValue<T> {
        HURCHALLA_FORCE_INLINE C() = default;
        // implicitly convert from canonical value C to montgomery value V
        HURCHALLA_FORCE_INLINE operator V() const
        {
            constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
            HPBC_PRECONDITION2(0 <= this->get() && this->get() < Rdiv2);
            return V(static_cast<SignedT>(this->get()));
        }
        HURCHALLA_FORCE_INLINE friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        HURCHALLA_FORCE_INLINE friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }
     protected:
        template <typename> friend class MontyHalfRange;
        template <template<class> class, template<class> class, typename>
          friend class MontyCommonBase;
        HURCHALLA_FORCE_INLINE explicit C(T a) : BaseMontgomeryValue<T>(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    struct FV : public V {
        HURCHALLA_FORCE_INLINE FV() = default;
     protected:
        template <typename> friend class MontyHalfRange;
        HURCHALLA_FORCE_INLINE explicit FV(SignedT a) : V(a) {}
    };
};


// Let the theoretical constant R = 1<<(ut_numeric_limits<T>::digits).
template <typename T>
class MontyHalfRange final :
                  public MontyCommonBase<MontyHalfRange, MontyHRValueTypes, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    static_assert(ut_numeric_limits<T>::digits >= 1, "");

    using S = typename MontyHRValueTypes<T>::SignedT;
    static_assert(ut_numeric_limits<S>::is_signed, "");
#if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<S>(-1) >> 1) == static_cast<S>(-1),
                          "Arithmetic right shift is required but unavailable");
#endif
    static_assert(static_cast<S>(-1) == ~(static_cast<S>(0)),
                                  "S must use two's complement representation");
    static_assert(static_cast<S>(static_cast<T>(static_cast<S>(-1))) ==
                  static_cast<S>(-1), "Casting a signed S value to unsigned and"
                               " back again must result in the original value");

    using BC = MontyCommonBase<::hurchalla::detail::MontyHalfRange,
                               ::hurchalla::detail::MontyHRValueTypes, T>;
    using BC::n_;
    using typename BC::V;
    using typename BC::C;
    using FV = typename MontyHRValueTypes<T>::FV;
 public:
    using MontyTag = TagMontyHalfrange;
    using uint_type = T;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;

    explicit MontyHalfRange(T modulus) : BC(modulus)
    {
        // MontyHalfRange requires  modulus < R/2
        constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(modulus < Rdiv2);
    }
    MontyHalfRange(const MontyHalfRange&) = delete;
    MontyHalfRange& operator=(const MontyHalfRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1)) - 1);
    }

    HURCHALLA_FORCE_INLINE V negate(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        S sx = x.get();
        S result = static_cast<S>(-sx);
        HPBC_POSTCONDITION2(isValid(V(result)));
        HPBC_POSTCONDITION2(getCanonicalValue(V(result)) ==
                                          getCanonicalValue(subtract(C(0), x)));
        return V(result);
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_INVARIANT2(n_ < Rdiv2);
        S result = halfrange_get_canonical<S>::call(x.get(),static_cast<S>(n_));
        HPBC_POSTCONDITION2(0 <= result && result < static_cast<S>(n_));
        return C(static_cast<T>(result));
    }

    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_INVARIANT2(n_ < Rdiv2);
        C cx = getCanonicalValue(x);
        T a = cx.get();
        HPBC_ASSERT2(a < n_);  // a has range [0, n_)

        HPBC_INVARIANT2(n_ % 2 == 1);
        T half_n_floor = static_cast<T>(n_ >> 1);     // == (n_ - 1) / 2

#if !defined(HURCHALLA_AVOID_CSELECT)
        T tmp = static_cast<T>(a - n_);
        // static_cast<S>(tmp) has range [-n, 0)
        HPBC_ASSERT2(-static_cast<S>(n_) <= static_cast<S>(tmp)
                     && static_cast<S>(tmp) < 0);
          // result = (half_n_floor < a) ? tmp : a
        T result = ::hurchalla::conditional_select((half_n_floor < a), tmp, a);
#else
        bool cond = (half_n_floor < a);
        T mask = static_cast<T>(-static_cast<T>(cond));
        T masked_n = static_cast<T>(mask & n_);
        T result = static_cast<T>(a - masked_n);
#endif
        // Assume a <= half_n_floor; then trivially we would have
        //    0 <= a <= (n_ - 1) / 2
        //    Given our assumption, we would have  0 <= result <= (n_ - 1) / 2
        // Assume a > half_n_floor; then we would have
        //    static_cast<S>(a - n_) > static_cast<S>(half_n_floor - n_)
        //    static_cast<S>(tmp) > static_cast<S>(((n_ - 1) / 2) - n_)
        //    static_cast<S>(tmp) >= static_cast<S>(((n_ - 1) / 2) - n_ + 1)
        //    static_cast<S>(tmp) >= -static_cast<S>((n_ - 1) / 2)
        //    From an assertion above, we already know static_cast<S>(tmp) < 0
        //    and so together we would have
        //    -static_cast<S>((n_ - 1) / 2) <= static_cast<S>(tmp) < 0.
        //    Therefore given our assumption, we would have
        //    -static_cast<S>((n_ - 1) / 2) <= static_cast<S>(result) < 0.
        // Since one of the two assumptions must be true, we have a possible
        // range for result of
        // -static_cast<S>((n_-1)/2) <= static_cast<S>(result) <= (n_-1)/2
        HPBC_POSTCONDITION2(
                     -static_cast<S>(half_n_floor) <= static_cast<S>(result)
                     && static_cast<S>(result) <= static_cast<S>(half_n_floor));
        return FV(static_cast<S>(result));
    }

    using BC::fmsub;
    // Multiplies two montgomery values x and y, and then subtracts the fusing-
    // value fv from the product.  Returns the resulting montgomery value.
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, FV fv, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_INVARIANT2(n_ % 2 == 1);
        // note that the constructor also established the invariant n < R/2.
        HPBC_PRECONDITION2(-(static_cast<S>((n_ - 1)/2)) <= fv.get()
                           && fv.get() <= static_cast<S>((n_ - 1)/2));
        namespace hc = ::hurchalla;
        T u_lo;
        S u_hi = hc::signed_multiply_to_hilo_product(u_lo, x.get(), y.get());

        // Performing the modular sub prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // subtraction.  See fmadd() in MontyCommonBase.h for a proof which
        // could be adapted to the subtraction in here.

        // By isValid(), x and y satisfy  -n <= x < n  and  -n <= y < n
        // Thus x*y <= n*n < n*R/2 == ((n-1)/2)*R + R/2
        // And since x*y == u_hi*R + u_lo, we know  u_hi <= ((n-1)/2).
        // Also,  x*y > -n*n > -n*R/2 == -((n+1)/2)*R + R/2.
        // And since x*y == u_hi*R + u_lo, we know  u_hi >= -((n+1)/2).
        // Putting this together, we know:  -((n+1)/2) <= u_hi <= ((n-1)/2).
        HPBC_ASSERT2(-(static_cast<S>((n_+1)/2)) <= u_hi
                     && u_hi <= static_cast<S>((n_-1)/2));
        // Since we have precondition that  -((n-1)/2) <= fv.get() <= (n-1)/2,
        // we know that
        // -((n+1)/2) - (n-1)/2 <= (u_hi - fv.get()) <= ((n-1)/2) - (-((n-1)/2))
        // Thus  -n <= (u_hi - fv.get()) <= n - 1 < n.
        // Also, since n < R/2, we have  -R/2 < (u_hi - fv.get()) < R/2,  which
        // means (u_hi - fv.get()) should never overflow type S.
        u_hi = static_cast<S>(u_hi - fv.get());
        HPBC_ASSERT2(-static_cast<S>(n_) <= u_hi && u_hi < static_cast<S>(n_));

        V vu_hi(u_hi);
        HPBC_ASSERT2(isValid(vu_hi));
        C cu_hi = getCanonicalValue(vu_hi);
        T tu_hi = cu_hi.get();
        HPBC_ASSERT2(tu_hi < n_);

        V result = montyREDC(tu_hi, u_lo, PTAG());

        HPBC_POSTCONDITION2(isValid(result));
        return result;
    }

    using BC::fmadd;
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, FV fv, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_INVARIANT2(n_ % 2 == 1);
        // note that the constructor also established the invariant n < R/2.
        HPBC_PRECONDITION2(-(static_cast<S>((n_ - 1)/2)) <= fv.get()
                           && fv.get() <= static_cast<S>((n_ - 1)/2));
        namespace hc = ::hurchalla;
        T u_lo;
        S u_hi = hc::signed_multiply_to_hilo_product(u_lo, x.get(), y.get());

        // Performing the modular add prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // addition.  See fmadd() in MontyCommonBase.h for a proof which could
        // be adapted to the modular addition in here.

        // See fmsub() above for why the following assert is true
        HPBC_ASSERT2(-(static_cast<S>((n_+1)/2)) <= u_hi
                     && u_hi <= static_cast<S>((n_-1)/2));
        // Since we have precondition that  -((n-1)/2) <= fv.get() <= (n-1)/2,
        // we know that
        // -((n+1)/2) + -((n-1)/2) <= (u_hi + fv.get()) <= ((n-1)/2) + (n-1)/2
        // Thus  -n <= (u_hi + fv.get()) <= n - 1 < n.
        // Also, since n < R/2, we have  -R/2 < (u_hi + fv.get()) < R/2,  which
        // means (u_hi + fv.get()) should never overflow type S.
        u_hi = static_cast<S>(u_hi + fv.get());
        HPBC_ASSERT2(-static_cast<S>(n_) <= u_hi && u_hi < static_cast<S>(n_));

        V vu_hi(u_hi);
        HPBC_ASSERT2(isValid(vu_hi));
        C cu_hi = getCanonicalValue(vu_hi);
        T tu_hi = cu_hi.get();
        HPBC_ASSERT2(tu_hi < n_);

        V result = montyREDC(tu_hi, u_lo, PTAG());

        HPBC_POSTCONDITION2(isValid(result));
        return result;
    }

    HURCHALLA_FORCE_INLINE V add(V x, C cy) const
    {
        HPBC_ASSERT2(isValid(x));        // we know  -n <= x.get() < n
        T tx = static_cast<T>(x.get());

#if !defined(HURCHALLA_AVOID_CSELECT)
        T tmpx = static_cast<T>(tx - n_);
        // Note if x.get() is negative, then tx >= n, since we know n < R/2.
        // This is assuming wrap-around conversion to T.  (ie two's complement?)
        // Likewise if tx >= n, then x.get() can't be >=0 because x.get() >= 0
        // would mean conversion to T would leave the value unchanged after
        // conversion: thus it would result in tx < n, which is obviously
        // impossible if tx >= n.  This means
        // that if tx >= n, then x.get() < 0.
        // And by contrapositive of the first item, tx < n implies x.get() >= 0
            // set tmpx = (tx>=n) ? tx : tmpx
        tmpx = ::hurchalla::conditional_select(tx >= n_, tx, tmpx);
#else
        // this code is functionally equivalent to the code in the #if above
        S sx = x.get();
        // Use arithmetic right shift of sign bit to create mask of all 1s or 0s
        T mask = static_cast<T>(sx >> ut_numeric_limits<S>::digits);
        // mask is all 1s if x < 0, and all 0s if x >= 0.
        T maskflip = static_cast<T>(~mask);
        // maskflip is all 1s if x >= 0, and all 0s if x < 0.
        T masked_n = static_cast<T>(maskflip & n_);
        // masked_n is n if x >= 0, and 0 if x < 0.
        T tmpx = static_cast<T>(tx - masked_n);
#endif
        // Assume tx >= n.  Then x.get() < 0, and the cselect will have set
        //    tmpx = tx = static_cast<T>(x.get()).
        //    And hence,  static_cast<S>(tmpx) == x.get() < 0.  We also know
        //    from isValid(x) that  -static_cast<S>(n) <= x.get().  Thus,
        //    -static_cast<S>(n) <= static_cast<S>(tmpx) < 0.
        // Assume tx < n.  Then x.get() >= 0, and the cselect will have left
        //    tmpx unchanged, so we know
        //    tmpx == tx - n == static_cast<T>(x.get()) - n,  and so
        //    static_cast<S>(tmpx) == x.get() - static_cast<S>(n).
        //    Since x.get() >= 0,  static_cast<S>(tmpx) >= -static_cast<S>(n).
        //    We know from isValid(x) that  x.get() < static_cast<S>(n),  and so
        //    static_cast<S>(tmpx) == x.get() - static_cast<S>(n) < n - n == 0.
        //    Putting this all together, we have
        //    -static_cast<S>(n) <= static_cast<S>(tmpx) < 0.
        // For both assumptions, we know we get
        // -static_cast<S>(n) <= static_cast<S>(tmpx) < 0.
        HPBC_ASSERT2(-static_cast<S>(n_) <= static_cast<S>(tmpx)
                     && static_cast<S>(tmpx) < 0);
        HPBC_ASSERT2(0 <= cy.get() && cy.get() < n_); // because cy is canonical
        S result = static_cast<S>(
                               static_cast<S>(tmpx) + static_cast<S>(cy.get()));
        // we can see result will satisfy  -n <= result < n - 1
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        return add(x, getCanonicalValue(y));
    }
    HURCHALLA_FORCE_INLINE C add(C cx, C cy) const
    {
        HPBC_PRECONDITION2(cx.get() < n_);
        HPBC_PRECONDITION2(cy.get() < n_);
        constexpr int max_bits_needed = ut_numeric_limits<T>::digits - 1;
        constexpr T limit = static_cast<T>(static_cast<T>(1)<<max_bits_needed);
        HPBC_ASSERT2(cx.get() < limit);
        HPBC_ASSERT2(cy.get() < limit);
        HPBC_ASSERT2(n_ < limit);
        S sx = static_cast<S>(cx.get());
        S sy = static_cast<S>(cy.get());
        S sn = static_cast<S>(n_);
        S modsum = ::hurchalla::modular_addition_prereduced_inputs(sx, sy, sn);
        HPBC_ASSERT2(modsum >= 0);
        T result = static_cast<T>(modsum);
        HPBC_POSTCONDITION2(result < n_);
        return C(result);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        C cx = getCanonicalValue(x);
        C cy = getCanonicalValue(y);
        S result = static_cast<S>(
                           static_cast<S>(cx.get()) - static_cast<S>(cy.get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract(V x, C cy) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(cy.get() < n_);
        C cx = getCanonicalValue(x);
        S result = static_cast<S>(
                           static_cast<S>(cx.get()) - static_cast<S>(cy.get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract(C cx, V y) const
    {
        HPBC_PRECONDITION2(cx.get() < n_);
        HPBC_PRECONDITION2(isValid(y));
        C cy = getCanonicalValue(y);
        S result = static_cast<S>(
                           static_cast<S>(cx.get()) - static_cast<S>(cy.get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE C subtract(C cx, C cy) const
    {
        HPBC_PRECONDITION2(cx.get() < n_);
        HPBC_PRECONDITION2(cy.get() < n_);
        constexpr int max_bits_needed = ut_numeric_limits<T>::digits - 1;
        constexpr T limit = static_cast<T>(static_cast<T>(1)<<max_bits_needed);
        HPBC_ASSERT2(cx.get() < limit);
        HPBC_ASSERT2(cy.get() < limit);
        HPBC_ASSERT2(n_ < limit);
        S sx = static_cast<S>(cx.get());
        S sy = static_cast<S>(cy.get());
        S sn = static_cast<S>(n_);
        namespace hc = ::hurchalla;
        S moddiff = hc::modular_subtraction_prereduced_inputs(sx, sy, sn);
        HPBC_ASSERT2(moddiff >= 0);
        T result = static_cast<T>(moddiff);
        HPBC_POSTCONDITION2(result < n_);
        return C(result);
    }

    template <typename J, typename K>
    HURCHALLA_FORCE_INLINE V unordered_subtract(J x, K y) const
    {
        return subtract(x, y);
    }

private:
    // functions called by the 'curiously recurring template pattern' base (BC).
    friend BC;

    // note: PTAG gets ignored here
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, PTAG) const
    {
        HPBC_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        namespace hc = ::hurchalla;
        bool isNegative;  // ignored
        T result = hc::REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
        resultIsZero = (result == 0);
        V v = V(static_cast<S>(result));
        HPBC_POSTCONDITION2(isValid(v));
        return v;
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(T u_hi, T u_lo, PTAG) const
    {
        bool resultIsZero;  // ignored
        return montyREDC(resultIsZero, u_hi, u_lo, PTAG());
    }

    // Let u be an arbitrary double-word value that is congruent (mod n_) to the
    // product of x and y, and that satisfies 0 <= u < n_*R.  Return the high
    // word of u, and write the low word of u to u_lo.
    HURCHALLA_FORCE_INLINE T multiplyToHiLo(T& u_lo, V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        S product_hi = ::hurchalla::signed_multiply_to_hilo_product(
                                                        u_lo, x.get(), y.get());
        // By isValid(), x and y both have range [-n, n).  Thus x*y has range
        // (-n*n, n*n], and due to class invariant n<R/2, x*y has range
        // (-n*R/2, n*R/2).  Since  product_hi*R + u_lo == x*y,
        // -n*R/2 < product_hi*R + u_lo < n*R/2.  Thus
        // product_hi*R < n*R/2 - u_lo <= n*R/2 < n*R/2 + R/2 == (n+1)*R/2
        // Since product_hi*R < ((n+1)/2)*R  and  n is odd, 
        // product_hi < (n + 1)/2 == (n-1)/2 + 1.  Thus product_hi <= (n-1)/2.
        // Likewise, since  product_hi*R + u_lo > -n*R/2,
        // product_hi*R > -n*R/2 - u_lo > -n*R/2 - R > -n*R/2 - 3*R/2,  and so
        // product_hi*R > -((n+3)/2)*R.  Therefore
        // product_hi > -((n+3)/2) == -((n+1)/2) - 1.  And thus
        // product_hi >= -(n+1)/2.  Putting it all together,
        // -(n+1)/2 <= product_hi <= (n-1)/2.  Also, due to invariant  n < R/2,
        // -R/2 < -(n+1)/2 <= product_hi <= (n-1)/2 < R/2.
#if 0
        T tmp = static_cast<T>(product_hi);
        T u_hi = static_cast<T>(tmp + n_);
        // We can easily deduce from the results in the above comments that
        // (n-1)/2 <= u_hi <= n + (n-1)/2 < R.

        // If (u_hi < n) is true, it indicates product_hi < 0.  Likewise,
        // (u_hi >= n) true would indicate product_hi >= 0.
           // Next line sets  uhi = (uhi>=n) ? tmp : uhi
        u_hi = ::hurchalla::conditional_select(u_hi >= n_, tmp, u_hi);
#else
        // this should be equivalent to the code above.  It's a slight hack
        // since product_hi doesn't actually represent a montgomery value (V),
        // but the spirit (and reality) of what getCanonicalValue() does is the
        // same as what we need, and using it lets us centralize optimal
        // handling of conditional-selects into getCanonicalValue().
        V v(product_hi);
        HPBC_ASSERT2(isValid(v));
        T u_hi = getCanonicalValue(v).get();
#endif
        HPBC_POSTCONDITION2(0 <= u_hi && u_hi < n_);
        return u_hi;
    }
    HURCHALLA_FORCE_INLINE T squareToHiLo(T& u_lo, V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        namespace hc = ::hurchalla;
        S tmp_hi = hc::signed_multiply_to_hilo_product(u_lo, x.get(), x.get());
        // The same logic as given in multiplyToHiLo shows that
        // -(n+1)/2 <= tmp_hi <= (n-1)/2.  But additionally, since the square
        // of an integer is always >= 0, we therefore know
        // 0 <= tmp_hi <= (n-1)/2.
        HPBC_ASSERT2(tmp_hi >= 0);
        T u_hi = static_cast<T>(tmp_hi);
        HPBC_POSTCONDITION2(0 <= u_hi && u_hi < n_);
        return u_hi;
    }
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        static_assert(!(ut_numeric_limits<T>::is_signed), "");
        constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_INVARIANT2(n_ < Rdiv2);
        return (-static_cast<S>(n_) <= x.get() && x.get() < static_cast<S>(n_));
    }
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        static_assert(std::is_same<S, decltype(x.get())>::value, "");
        return (0 <= x.get() && x.get() < static_cast<S>(n_));
    }
    // get a Natural number (i.e. number >= 0) congruent to x (mod n)
    HURCHALLA_FORCE_INLINE T getNaturalEquivalence(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        T result = static_cast<T>(static_cast<T>(x.get()) + n_);
        // Since the precondition of isValid(x) required
        // -(static_cast<S>(n)) <= x.get() < n,  the sum for 'result' will have
        // carried (wrapped-around) if  x.get() < 0.
        HPBC_POSTCONDITION2(0 <= result && result < static_cast<T>(2*n_));
        return result;
    }
};


}} // end namespace

#endif
