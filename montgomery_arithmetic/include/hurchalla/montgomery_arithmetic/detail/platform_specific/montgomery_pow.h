// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/Unroll.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <array>
#include <cstddef>
#include <type_traits>

namespace hurchalla { namespace detail {


// This file is intended to implement the class MontgomeryForm's member
// functions pow() and array_pow().

// MF should be a MontgomeryForm type
template<class MF>
struct montgomery_pow {
  using T = typename MF::IntegerType;
  using V = typename MF::MontgomeryValue;
  static_assert(ut_numeric_limits<T>::is_integer, "");

  static HURCHALLA_FORCE_INLINE V scalarpow(const MF& mf, V base, T exponent)
  {
    HPBC_PRECONDITION(exponent >= 0);
    // This is an optimized version of Algorithm 14.76, from
    // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
    // See also: hurchalla/modular_arithmetic/detail/impl_modular_pow.h
    V result;
    if (exponent & static_cast<T>(1))
        result = base;
    else
        result = mf.getUnityValue();
    while (exponent > static_cast<T>(1)) {
        exponent = static_cast<T>(exponent >> static_cast<T>(1));
        base = mf.template multiply<LowuopsTag>(base, base);
        // The multiply above is a loop carried dependency.  Thus, a second loop
        // carried dependency with the same length can be essentially free due
        // to instruction level parallelism, so long as it does not introduce
        // any branch mispredictions.
        // So we will always compute the second multiply, instead of
        // conditionally computing it, and we will encourage the compiler to use
        // a (branchless) conditional move instruction.
        // We use lowlatencyTag below because the 'result' loop carried
        // dependency depends upon both multiply and a conditional move, whereas
        // 'base' above depends only on multiply and thus is tagged for lowuops
        // since it is less likely to be a latency bottleneck.
        V tmp = mf.template multiply<LowlatencyTag>(result, base);
#if 0
        // ternary op generally compiles to conditional moves.
        result = (exponent & static_cast<T>(1)) ? tmp : result;
#else
        HURCHALLA_CMOV(exponent & static_cast<T>(1), result, tmp);
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
            tmp[i] = mf.template multiply<LowuopsTag>(result[i], bases[i]);
        });
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
#if 0
            //the ternary operator usually/hopefully results in our desired cmov
            result[i] = (exponent & static_cast<T>(1)) ? tmp[i] : result[i];
#else
            HURCHALLA_CMOV(exponent & static_cast<T>(1), result[i], tmp[i]);
#endif
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
        using U = decltype(result[0].value);
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
            V tmp = mf.template multiply<LowuopsTag>(result[i], bases[i]);
            result[i].value = (mask & tmp.value) | (maskflip & result[i].value);
        });
    }
    return result;
  }
};


// MF should be a MontgomeryForm type
template<class MF>
struct DefaultMontArrayPow {
    using T = typename MF::IntegerType;
    using V = typename MF::MontgomeryValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    // Below is a catch-all template version of pow(), followed by overloads of
    // pow() for different array sizes - by C++ rules, overloads take precedence
    // over templates whenever argument matching succeeds.

    // Template version (with two std::enable_if forms):
    // Having a std::enable_if cutoff at NUM_BASES < 10 is a bit arbitrary,
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
    typename std::enable_if<(NUM_BASES < 10), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, const std::array<V, NUM_BASES>& bases, T exponent)
    {
        static_assert(NUM_BASES > 0, "");
        // conditional branching seems to typically work best for large-ish
        // array sizes.  And so long as the array isn't huge, unrolling helps.
        return montgomery_pow<MF>::arraypow_cond_branch_unrolled(
                                                           mf, bases, exponent);
    }
    template <std::size_t NUM_BASES>
    static HURCHALLA_FORCE_INLINE
    typename std::enable_if<(NUM_BASES >= 10), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, const std::array<V, NUM_BASES>& bases, T exponent)
    {
        static_assert(NUM_BASES > 0, "");
        // When NUM_BASES gets huge we no longer want to force-unroll loops, but
        // other than that we'll do the same as arraypow_cond_branch_unrolled().
        return montgomery_pow<MF>::arraypow_cond_branch(mf, bases, exponent);
    }

    // delegate a size 1 array call to the scalar (non-array) montgomery_pow
    static HURCHALLA_FORCE_INLINE
    std::array<V,1> pow(const MF& mf, const std::array<V,1>& bases, T exponent)
    {
        std::array<V,1> result;
        result[0] = montgomery_pow<MF>::scalarpow(mf, bases[0], exponent);
        return result;
    }

    static HURCHALLA_FORCE_INLINE
    std::array<V,2> pow(const MF& mf, const std::array<V,2>& bases, T exponent)
    {
        return montgomery_pow<MF>::arraypow_cmov(mf, bases, exponent);
    }
};



// Primary template.  Forwards all work to DefaultMontArrayPow.
template<class MontyType, class MF, class Enable = void>
struct montgomery_array_pow {
    using V = typename MF::MontgomeryValue;
    template <std::size_t NUM_BASES>
    static HURCHALLA_FORCE_INLINE
    std::array<V, NUM_BASES> pow(const MF& mf,
                                 const std::array<V, NUM_BASES>& bases,
                                 typename MF::IntegerType exponent)
    {
        return DefaultMontArrayPow<MF>::pow(mf, bases, exponent);
    }
};

#if defined(HURCHALLA_TARGET_ISA_X86_64)
// Partial specialization using QuarterrangeTag (x86-64 only).
// I haven't measured performance for any ISAs except x86-64, and so all other
// ISAs default to the primary template above.
//
// We use the awkward "enable_if" syntax in this template because we want to
// match MontyType::MontyTag to QuarterrangeTag - but the client code can't
// provide MontyType::MontyTag because sometimes it might not exist and the
// client would get a compile error if so.  But in the template parameters below
// if MontyType::MontyTag doesn't exist that's ok, because due to SFINAE the
// compiler will just give up on attempting the match (it will instead match to
// the primary template above).
template<class MontyType, class MF>
struct montgomery_array_pow<MontyType, MF, typename std::enable_if<
              std::is_same<typename MontyType::MontyTag, QuarterrangeTag>::value
          >::type> {
    using T = typename MF::IntegerType;
    using V = typename MF::MontgomeryValue;

    template <std::size_t NUM_BASES>
    static HURCHALLA_FORCE_INLINE
    std::array<V, NUM_BASES> pow(const MF& mf,
                                 const std::array<V, NUM_BASES>& bases,
                                 T exponent)
    {
        return DefaultMontArrayPow<MF>::pow(mf, bases, exponent);
    }

    // For x86-64 Montgomery QuarterrangeTag and 2 or 3 bases and compiling with
    // clang, arraypow_cmov() always had better measured performance than any of
    // the alternative functions.  For gcc and icc compilers, arraypow_masked()
    // resulted in the best measured performance.  But in general we expect
    // conditional moves (cmovs) to perform better than masks, so we default for
    // any unmeasured compiler to using arraypow_cmov().  I haven't measured
    // perf on MSVC, for example.
#  if !defined(__GNUC__) || defined(__clang__)
    static HURCHALLA_FORCE_INLINE
    std::array<V,2> pow(const MF& mf, const std::array<V,2>& bases, T exponent)
    {
        return montgomery_pow<MF>::arraypow_cmov(mf, bases, exponent);
    }
    static HURCHALLA_FORCE_INLINE
    std::array<V,3> pow(const MF& mf, const std::array<V,3>& bases, T exponent)
    {
        return montgomery_pow<MF>::arraypow_cmov(mf, bases, exponent);
    }

#  else
    static HURCHALLA_FORCE_INLINE
    std::array<V,2> pow(const MF& mf, const std::array<V,2>& bases, T exponent)
    {
        return montgomery_pow<MF>::arraypow_masked(mf, bases, exponent);
    }
    static HURCHALLA_FORCE_INLINE
    std::array<V,3> pow(const MF& mf, const std::array<V,3>& bases, T exponent)
    {
        return montgomery_pow<MF>::arraypow_masked(mf, bases, exponent);
    }
#  endif
};
#endif


}} // end namespace

#endif
