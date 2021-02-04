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


// This class is intended solely for internal use by this file
// MF should be a MontgomeryForm type
template<class MF>
struct MontPowImpl {
  using T = typename MF::T_type;
  using V = typename MF::MontgomeryValue;
  static_assert(ut_numeric_limits<T>::is_integer, "");

  static HURCHALLA_FORCE_INLINE V pow(const MF& mf, V base, T exponent)
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
  // These array versions have a performance advantage due to instruction level
  // parallelism, compared to the non-array montgomery_pow function.
  // They use the same algorithm as the non-array montgomery_pow().
  // These array version functions should all be equivalent to one another,
  // aside from their differences in performance.
  // --------
  template <std::size_t NUM_BASES>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_cond_branch(const MF& mf, std::array<V, NUM_BASES> bases, T exponent)
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
        exponent = static_cast<T>(exponent >> static_cast<T>(1));
        T lowbit = (exponent & static_cast<T>(1));
        using U = decltype(result[0].get());
        static_assert(ut_numeric_limits<U>::is_integer, "");
        static_assert(!(ut_numeric_limits<U>::is_signed), "");
        U mask = static_cast<U>(0) - static_cast<U>(lowbit);
        U maskflip = static_cast<U>(lowbit) - static_cast<U>(1);
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template multiply<LowuopsTag>(bases[i], bases[i]);
// TODO: use LowuopsTag?
            V tmp = mf.template multiply<LowlatencyTag>(result[i], bases[i]);
            result[i]= V((mask & (tmp.get())) | (maskflip & (result[i].get())));
        });
    }
    return result;
  }

  // This version of the function version is intended for very large NUM_BASES.
  // For example, on haswell NUM_BASES of 96 performs quite efficiently,
  // especially with compiler flags -march=haswell and -mtune=haswell.
  // Note that being able to find a huge number of bases that you want to use
  // for modular exponentiation all with the same power and modulus, is likely a
  // special case.  Even just having more than one single base, so that you can
  // use an arraypow function rather than pow, may be unusual.
  template <std::size_t NUM_BASES>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_huge(const MF& mf, std::array<V, NUM_BASES> bases, T exponent)
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
};


// delegation class
// MF should be a MontgomeryForm type
template<class MF>
struct MontPow {
    using T = typename MF::T_type;
    using V = typename MF::MontgomeryValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");

    // The catch-all template version.
    // Having the std::enable_if cutoff at NUM_BASES < 32 is a bit arbitrary,
    // since the only way to know the best cutoff for any given machine is to
    // make performance measurements.  However, it may be roughly okay in
    // practice since most of the time NUM_BASES is likely to be < 10, and
    // otherwise (probably much less common), it's likely > 50, in the "huge"
    // category where arraypow_huge() will likely perform best.  Even if
    // arraypow_huge() doesn't perform best in microbenchmarks, its much smaller
    // function size compared to arraypow_cond_branch (and consequent reduced
    // i-cache use) might provide more of a benefit to the whole program's
    // performance as NUM_BASES starts to get "huge".
    template <std::size_t NUM_BASES>
    static HURCHALLA_FORCE_INLINE
    typename std::enable_if<(NUM_BASES < 32), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, std::array<V, NUM_BASES>& bases, T exponent)
    {
        static_assert(NUM_BASES > 0, "");
        // conditional branching seems to typically work best for large-ish
        // array sizes.
        return MontPowImpl<MF>::arraypow_cond_branch(mf, bases, exponent);
    }

    template <std::size_t NUM_BASES>
    static HURCHALLA_FORCE_INLINE
    typename std::enable_if<(NUM_BASES >= 32), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, std::array<V, NUM_BASES>& bases, T exponent)
    {
        static_assert(NUM_BASES > 0, "");
        // When NUM_BASES gets huge we don't want to force-unroll loops, but
        // other than that we do the same as arraypow_cond_branch
        return MontPowImpl<MF>::arraypow_huge(mf, bases, exponent);
    }

    // delegate a size 1 array call to the non-array montgomery_pow
    static HURCHALLA_FORCE_INLINE
    std::array<V,1> pow(const MF& mf, std::array<V,1>& bases, T exponent)
    {
        HPBC_PRECONDITION(exponent >= 0);
        std::array<V,1> result;
        result[0] = MontPowImpl<MF>::pow(mf, bases[0], exponent);
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

#if defined(HURCHALLA_TARGET_ISA_X86_64)
// I've only measured x86-64 so far - I measured this section to be a perf
// improvement over the catch-all template version.  In general I'd expect the
// template version to work best for an array sized 3 or larger though, so I've
// only enabled this section for the (measured) x86-64 ISA.
#  if defined(__GNUC__) && !defined(__clang__)
    // This section will be enabled for both gcc and icc, but not clang.
    // x86-64 gcc seems to do better with masking for array size of 3.  But in
    // general we expect conditional moves (cmovs) to perform better than masks.
    // x86-64 icc does much better with the masking version (at size 3) than
    // with the cmov version.
    static HURCHALLA_FORCE_INLINE
    std::array<V,3> pow(const MF& mf, std::array<V,3>& bases, T exponent)
    {
        return MontPowImpl<MF>::arraypow_masked(mf, bases, exponent);
    }
#  else
    static HURCHALLA_FORCE_INLINE
    std::array<V,3> pow(const MF& mf, std::array<V,3>& bases, T exponent)
    {
        return MontPowImpl<MF>::arraypow_cmov(mf, bases, exponent);
    }
#  endif
#endif
};


template <class MF>
HURCHALLA_FORCE_INLINE typename MF::MontgomeryValue
montgomery_pow(const MF& mf, typename MF::MontgomeryValue base,
              typename MF::T_type exponent)
{
    return MontPowImpl<MF>::pow(mf, base, exponent);
}

template <class MF, std::size_t NUM_BASES>
HURCHALLA_FORCE_INLINE std::array<typename MF::MontgomeryValue, NUM_BASES>
montgomery_pow(const MF& mf,
               std::array<typename MF::MontgomeryValue, NUM_BASES>& bases,
               typename MF::T_type exponent)
{
    return MontPow<MF>::pow(mf, bases, exponent);
}


}} // end namespace

#endif
