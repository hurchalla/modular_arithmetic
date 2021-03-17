// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/Unroll.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <array>
#include <cstddef>

namespace hurchalla { namespace detail {


// The montgomery_pow() function at bottom of this file implements the class
// MontgomeryForm's member function pow().

// This MontPowImpl struct is intended solely for internal use by this file.
// MF should be a MontgomeryForm type
template<class MF>
struct MontPowImpl {
  using T = typename MF::T_type;
  using V = typename MF::MontgomeryValue;
  static_assert(ut_numeric_limits<T>::is_integer, "");

  static HURCHALLA_FORCE_INLINE V scalarpow(const MF& mf, V base, T exponent)
  {
    HPBC_PRECONDITION(exponent >= 0);
    // This is an optimized version of Algorithm 14.76, from
    // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
    // See also: hurchalla/modular_arithmetic/detail/impl_modular_pow.h
    V result = (exponent & static_cast<T>(1)) ? base : mf.getUnityValue();
    while (exponent > static_cast<T>(1)) {
        exponent = static_cast<T>(exponent >> static_cast<T>(1));
        base = mf.template multiply<LowuopsTag>(base, base);
        // The multiply above is a loop carried dependency.  Thus, a second loop
        // carried dependency with the same length can be essentially free due
        // to instruction level parallelism, so long as it does not introduce
        // ny branch mispredictions.
        // So we will always compute the second multiply, instead of
        // conditionally computing it, and we will encourage the compiler to use
        // a (branchless) conditional move instruction.
        // We use lowlatencyTag below because the 'result' loop carried
        // dependency depends upon both multiply and a conditional move, whereas
        // 'base' above depends only on multiply and thus is tagged for lowuops
        // since it is less likely to be a latency bottleneck.
        V tmp = mf.template multiply<LowlatencyTag>(result, base);
#if 1
  // ternary op generally compiles to conditional moves.  On x64 gcc and clang
  // performance was significantly better with this line than with the masking
  // method below.  I haven't measured icc or msvc.
        result = (exponent & static_cast<T>(1)) ? tmp : result;
#else
        T lowbit = (exponent & static_cast<T>(1));
        using U = decltype(tmp.get());
        static_assert(ut_numeric_limits<U>::is_integer, "");
        static_assert(!(ut_numeric_limits<U>::is_signed), "");
        U mask = static_cast<U>(0) - static_cast<U>(lowbit);
        U mask_flipped = static_cast<U>(lowbit) - static_cast<U>(1);
        result = V((mask & (tmp.get())) | (mask_flipped & (result.get())));
#endif
    }
    return result;
  }

  // --------
  // The array versions below have a performance advantage due to instruction
  // level parallelism, compared to the non-array montgomery_pow function.
  // They use the same algorithm as the non-array montgomery_pow().  They
  // require an array of bases, of which each element is modularly exponentiated
  // to the same power and returned as an array.  At least one application is
  // miller-rabin primality testing, since it needs to raise multiple bases to
  // the same power.
  // These array version functions should all be equivalent to one another,
  // aside from their differences in performance.
  // --------

  // This first arraypow function version is the most obvious implementation and
  // results in the smallest code size.  It works well when you wish to limit
  // the code size and/or the i-cache usage, but it also is competitive at
  // providing the best performance when NUM_BASES gets huge (roughly values of
  // NUM_BASES > 50, though note that such situations are probably unusual in
  // practice, given that even having NUM_BASES > 1 could be considered a
  // special case for pow).
  template <std::size_t NUM_BASES>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_cond_branch(const MF& mf, std::array<V, NUM_BASES> bases, T exponent)
  {
    HPBC_PRECONDITION(exponent >= 0);
    std::array<V, NUM_BASES> result;
    if (exponent & static_cast<T>(1)) {
        for (std::size_t i=0; i<NUM_BASES; ++i)
            result[i] = bases[i];
    } else {
        V unity = mf.getUnityValue();
        for (std::size_t i=0; i<NUM_BASES; ++i)
            result[i] = unity;
    }
    while (exponent > static_cast<T>(1)) {
        exponent = static_cast<T>(exponent >> static_cast<T>(1));
        for (std::size_t i=0; i<NUM_BASES; ++i)
            bases[i] = mf.template multiply<LowuopsTag>(bases[i], bases[i]);
        if (exponent & static_cast<T>(1)) {
            for (std::size_t i=0; i<NUM_BASES; ++i)
                result[i]= mf.template multiply<LowuopsTag>(result[i],bases[i]);
        }
    }
    return result;
  }

  template <std::size_t NUM_BASES>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_cond_branch_unrolled(const MF& mf, std::array<V, NUM_BASES> bases,
                                                                     T exponent)
  {
    HPBC_PRECONDITION(exponent >= 0);
    std::array<V, NUM_BASES> result;
    if (exponent & static_cast<T>(1)) {
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            result[i] = bases[i];
        });
    } else {
        V unity = mf.getUnityValue();
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            result[i] = unity;
        });
    }
    while (exponent > static_cast<T>(1)) {
        exponent = static_cast<T>(exponent >> static_cast<T>(1));
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template multiply<LowuopsTag>(bases[i], bases[i]);
        });
        if (exponent & static_cast<T>(1)) {
            Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
                result[i]= mf.template multiply<LowuopsTag>(result[i],bases[i]);
            });
        }
    }
    return result;
  }

  template <std::size_t NUM_BASES>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_cmov(const MF& mf, std::array<V, NUM_BASES> bases, T exponent)
  {
    HPBC_PRECONDITION(exponent >= 0);
    std::array<V, NUM_BASES> result;
// TODO: use cmov/masking on this?
    if (exponent & static_cast<T>(1)) {
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            result[i] = bases[i];
        });
    } else {
        V unity = mf.getUnityValue();
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            result[i] = unity;
        });
    }
    while (exponent > static_cast<T>(1)) {
        exponent = static_cast<T>(exponent >> static_cast<T>(1));
        std::array<V, NUM_BASES> tmp;
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template multiply<LowuopsTag>(bases[i], bases[i]);
// TODO: use LowuopsTag?
            tmp[i] = mf.template multiply<LowlatencyTag>(result[i], bases[i]);
        });
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            //the ternary operator usually/hopefully results in our desired cmov
            result[i] = (exponent & static_cast<T>(1)) ? tmp[i] : result[i];
        });
    }
    return result;
  }

  template <std::size_t NUM_BASES>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_masked(const MF& mf, std::array<V, NUM_BASES> bases, T exponent)
  {
    HPBC_PRECONDITION(exponent >= 0);
    std::array<V, NUM_BASES> result;
// TODO: use cmov/masking on this?
    if (exponent & static_cast<T>(1)) {
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            result[i] = bases[i];
        });
    } else {
        V unity = mf.getUnityValue();
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            result[i] = unity;
        });
    }
    while (exponent > static_cast<T>(1)) {
        using U = decltype(result[0].get());
        static_assert(ut_numeric_limits<U>::is_integer, "");
        static_assert(!(ut_numeric_limits<U>::is_signed), "");
        static_assert(ut_numeric_limits<U>::digits >=
                      ut_numeric_limits<T>::digits, "");
        exponent = static_cast<T>(exponent >> 1u);
        U lowbit = static_cast<U>(static_cast<U>(exponent) & static_cast<U>(1));
        U mask = static_cast<U>(static_cast<U>(0) - lowbit);
        U maskflip = static_cast<U>(lowbit - static_cast<U>(1));
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template multiply<LowuopsTag>(bases[i], bases[i]);
// TODO: use LowuopsTag?
            V tmp = mf.template multiply<LowlatencyTag>(result[i], bases[i]);
            result[i]= V((mask & (tmp.get())) | (maskflip & (result[i].get())));
        });
    }
    return result;
  }
};


// MF should be a MontgomeryForm type
template<class MF>
struct MontArrayPow {
    using T = typename MF::T_type;
    using V = typename MF::MontgomeryValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    // Below is a catch-all template version of pow(), followed by overloads of
    // pow() for different array sizes - by C++ rules, overloads take precedence
    // over templates whenever argument matching succeeds.

    // Template version (with two std::enable_if forms):
    // Having a std::enable_if cutoff at NUM_BASES < 20 is a bit arbitrary,
    // since the only way to know the best cutoff for any given machine is to
    // make performance measurements.  However, it may be roughly okay in
    // practice since most of the time NUM_BASES is likely to be < 10, and
    // otherwise (probably much less common), it's likely > 50, in the "huge"
    // category where arraypow_cond_branch() will likely perform best.  Even if
    // arraypow_cond_branch() doesn't perform best in microbenchmarks, its much
    // smaller function size compared to arraypow_cond_branch_unrolled (and
    // consequent reduced i-cache use) might provide more of a benefit to the
    // whole program's performance as NUM_BASES starts to get "huge".
    // An additional reason that we would prefer to use a conservatively low
    // cutoff is that clang in debug mode uses up a huge amount of memory
    // (enough to crash a system while testing sometimes) during compilation of
    // the Unroll class when the NUM_BASES parameter is large.  Since the plain
    // arraypow_cond_branch() doesn't use Unroll, we avoid the problem by using
    // that function when NUM_BASES starts to get large.
    template <std::size_t NUM_BASES>
    static HURCHALLA_FORCE_INLINE
    typename std::enable_if<(NUM_BASES < 20), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, std::array<V, NUM_BASES>& bases, T exponent)
    {
        static_assert(NUM_BASES > 0, "");
        // conditional branching seems to typically work best for large-ish
        // array sizes.  And so long as the array isn't huge, unrolling helps.
        return MontPowImpl<MF>::arraypow_cond_branch_unrolled(
                                                           mf, bases, exponent);
    }
    template <std::size_t NUM_BASES>
    static HURCHALLA_FORCE_INLINE
    typename std::enable_if<(NUM_BASES >= 20), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, std::array<V, NUM_BASES>& bases, T exponent)
    {
        static_assert(NUM_BASES > 0, "");
        // When NUM_BASES gets huge we no longer want to force-unroll loops, but
        // other than that we'll do the same as arraypow_cond_branch_unrolled().
        return MontPowImpl<MF>::arraypow_cond_branch(mf, bases, exponent);
    }

    // delegate a size 1 array call to the scalar (non-array) montgomery_pow
    static HURCHALLA_FORCE_INLINE
    std::array<V,1> pow(const MF& mf, std::array<V,1>& bases, T exponent)
    {
        HPBC_PRECONDITION(exponent >= 0);
        std::array<V,1> result;
        result[0] = MontPowImpl<MF>::scalarpow(mf, bases[0], exponent);
        return result;
    }

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && \
    defined(HURCHALLA_TARGET_ISA_X86_64)
// x86-64 gcc seems to do better with masking for array size of 2.  But in
// general we expect conditional moves (cmovs) to perform better than masks.
// x86-64 icc has no difference between mask and cmov, at -O2 and -O3.  However
// at -Os (and -O1) icc takes ~4x as long as it should with cmov and ~3x as long
// as it should with mask.
//TODO: file a bug report regarding icc perf at -Os?  FYI using the compile flag
// -inline-forceinline gets it to ~1.5x the runtime of -O2 (which is acceptable)
// so this strongly suggests it's due to inlining not happening, but my VTune
// and godbolt experiments haven't yet shown inlining to be skipped at -Os when
// a function has __attribute__((always_inline)).  And I double checked to be
// certain that all involved functions had this (via HURCHALLA_FORCE_INLINE).
    static HURCHALLA_FORCE_INLINE
    std::array<V,2> pow(const MF& mf, std::array<V,2>& bases, T exponent)
    {
        return MontPowImpl<MF>::arraypow_masked(mf, bases, exponent);
    }
#else
    static HURCHALLA_FORCE_INLINE
    std::array<V,2> pow(const MF& mf, std::array<V,2>& bases, T exponent)
    {
        return MontPowImpl<MF>::arraypow_cmov(mf, bases, exponent);
    }
#endif

// I've only measured x86-64 so far - I measured the next section to be a perf
// improvement over the catch-all template version.  In general I'd expect the
// template version to work best for an array sized 3 or larger though, so I've
// only enabled the section for the (measured) x86-64 ISA.
#if defined(HURCHALLA_TARGET_ISA_X86_64)
#  if !defined(__GNUC__) || defined(__clang__)
    // In general we expect conditional moves (cmovs) to perform better than
    // masks, so we default to using arraypow_cmov rather than arraypow_masked.
    // Gcc and icc are exceptions; x86-64 icc does much better for array size
    // 3 with the masking version than with the cmov version, and x86-64 gcc
    // does slightly better with masking than cmov for size 3.
    static HURCHALLA_FORCE_INLINE
    std::array<V,3> pow(const MF& mf, std::array<V,3>& bases, T exponent)
    {
        return MontPowImpl<MF>::arraypow_cmov(mf, bases, exponent);
    }
#  else
    static HURCHALLA_FORCE_INLINE
    std::array<V,3> pow(const MF& mf, std::array<V,3>& bases, T exponent)
    {
        return MontPowImpl<MF>::arraypow_masked(mf, bases, exponent);
    }
#  endif
#endif
};



template <class MF>
HURCHALLA_FORCE_INLINE typename MF::MontgomeryValue
montgomery_pow(const MF& mf, typename MF::MontgomeryValue base,
              typename MF::T_type exponent)
{
    return MontPowImpl<MF>::scalarpow(mf, base, exponent);
}

template <class MF, std::size_t NUM_BASES>
HURCHALLA_FORCE_INLINE std::array<typename MF::MontgomeryValue, NUM_BASES>
montgomery_pow(const MF& mf,
               std::array<typename MF::MontgomeryValue, NUM_BASES>& bases,
               typename MF::T_type exponent)
{
    return MontArrayPow<MF>::pow(mf, bases, exponent);
}


}} // end namespace

#endif
