// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/Unroll.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

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
    HPBC_CLOCKWORK_PRECONDITION(exponent >= 0);
    // This is an optimized version of Algorithm 14.76, from
    // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
    // See also: hurchalla/modular_arithmetic/detail/impl_modular_pow.h
    V mont_one = mf.getUnityValue();
    V result = mont_one;
    result.cmov((exponent & static_cast<T>(1)), base);
    while (exponent > static_cast<T>(1)) {
        exponent = static_cast<T>(exponent >> static_cast<T>(1));
        base = mf.template square<LowlatencyTag>(base);
#if 0
        // The square above is a loop carried dependency.  Thus, a second loop
        // carried dependency with the same length can be essentially free due
        // to instruction level parallelism, so long as it does not introduce
        // any branch mispredictions.
        // So we will always compute the following multiply, instead of
        // conditionally computing it, and we will encourage the compiler to use
        // a (branchless) conditional move instruction.
        // We use lowlatencyTag below because the 'result' loop carried
        // dependency depends upon both multiply and a conditional move, whereas
        // 'base' above depends only on square (multiply) and thus is tagged for
        // lowuops since it is less likely to be a latency bottleneck.

        V tmp = mf.template multiply<LowlatencyTag>(result, base);
        result.cmov(exponent & static_cast<T>(1), tmp);
#else
        // This section shortens the dependency chain for 'result'.  Since
        // 'result' had a longer dependency chain above than 'base' did, this
        // in theory should run faster.  An in practice it has consistently
        // benchmarked better - usually ~5% faster, and never slower than above.
        V tmp = mont_one;
        tmp.cmov(exponent & static_cast<T>(1), base);
        result = mf.template multiply<LowlatencyTag>(result, tmp);
#endif
    }
    return result;
  }

  // --------
  // The array versions below have a performance advantage due to instruction
  // level parallelism, compared to the non-array montgomery_pow function.
  // They use the same algorithm as the non-array montgomery_pow().  They
  // require an array of bases, of which each element is modularly exponentiated
  // to the same power and returned as an array.  An example application where
  // array pow can be useful is miller-rabin primality testing, since it needs
  // to raise multiple bases to the same power.
  // These array version functions should all be equivalent to one another,
  // aside from their differences in performance.
  // --------

  // This first arraypow function version is the most obvious implementation and
  // results in the smallest code size.  It works well when you wish to limit
  // the code size and/or the i-cache usage, but it also is competitive at
  // providing the best performance when NUM_BASES gets huge (roughly values of
  // NUM_BASES > 50, though note that such situations are probably unusual in
  // practice, given that even having NUM_BASES > 1 could be considered a
  // special case for pow).  The best reason to use this version is likely to be
  // to avoid code bloat - when it has very little performance justification.
  template <std::size_t NUM_BASES, class PTAG>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_cond_branch(const MF& mf, std::array<V, NUM_BASES> bases, T exponent, PTAG)
  {
    HPBC_CLOCKWORK_PRECONDITION(exponent >= 0);
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
            bases[i] = mf.template square<PTAG>(bases[i]);
        if (exponent & static_cast<T>(1)) {
            for (std::size_t i=0; i<NUM_BASES; ++i)
                result[i]= mf.template multiply<PTAG>(result[i],bases[i]);
        }
    }
    return result;
  }

  template <std::size_t NUM_BASES, class PTAG>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_cond_branch_unrolled(const MF& mf, std::array<V, NUM_BASES> bases,
                                                               T exponent, PTAG)
  {
    HPBC_CLOCKWORK_PRECONDITION(exponent >= 0);
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
            bases[i] = mf.template square<PTAG>(bases[i]);
        });
        if (exponent & static_cast<T>(1)) {
            Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
                result[i] = mf.template multiply<PTAG>(result[i], bases[i]);
            });
        }
    }
    return result;
  }

  template <std::size_t NUM_BASES, class PTAG>
  static HURCHALLA_FORCE_INLINE std::array<V, NUM_BASES>
  arraypow_cmov(const MF& mf, std::array<V, NUM_BASES> bases, T exponent, PTAG)
  {
    HPBC_CLOCKWORK_PRECONDITION(exponent >= 0);
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
#if 0
        std::array<V, NUM_BASES> tmp;

        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template square<PTAG>(bases[i]);
            tmp[i] = mf.template multiply<PTAG>(result[i], bases[i]);
        });
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            result[i].cmov(exponent & static_cast<T>(1), tmp[i]);
        });
#else
        // see scalarpow)() comments for why this #else section might be
        // preferable to the #if alternative above.  There's probably little
        // difference at larger NUM_BASES though, where total uops is the
        // bottleneck rather than dependency chain length.
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template square<PTAG>(bases[i]);
        });

        V mont_one = mf.getUnityValue();
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            V tmp = mont_one;
            tmp.cmov(exponent & static_cast<T>(1), bases[i]);
            result[i] = mf.template multiply<PTAG>(result[i], tmp);
        });
#endif
    }
    return result;
  }

  template <std::size_t NUM_BASES, class PTAG>
  static HURCHALLA_FORCE_INLINE std::array<V,NUM_BASES>
  arraypow_masked(const MF& mf, std::array<V,NUM_BASES> bases, T exponent, PTAG)
  {
    HPBC_CLOCKWORK_PRECONDITION(exponent >= 0);
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
        exponent = static_cast<T>(exponent >> 1u);
#if 0
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template square<PTAG>(bases[i]);
            V tmp = mf.template multiply<PTAG>(result[i], bases[i]);
            result[i].template
                      cmov<CSelectMaskedTag>(exponent & static_cast<T>(1), tmp);
        });
#else
        // the comments in arraypow_cmov() apply equally in this section
        V mont_one = mf.getUnityValue();
        Unroll<NUM_BASES>::call([&](std::size_t i) HURCHALLA_INLINE_LAMBDA {
            bases[i] = mf.template square<PTAG>(bases[i]);
            V tmp = mont_one;
            tmp.template
                 cmov<CSelectMaskedTag>(exponent & static_cast<T>(1), bases[i]);
            result[i] = mf.template multiply<PTAG>(result[i], tmp);
        });
#endif
    }
    return result;
  }
};



// Primary template.
// This struct and its specializations are intended for NUM_BASES <= 3.
// There are partial specializations below for TagMontyHalfrange and
// TagMontyQuarterrange.  This primary template is mostly used for handling
// TagMontyFullrange.  (it would also be used for TagMontyWrappedmath)
template<class MontyTag, class MF>
struct montgomery_array_pow_small {
    using V = typename MF::MontgomeryValue;

    template <std::size_t NUM_BASES, class PTAG>
    static HURCHALLA_FORCE_INLINE
    std::array<V, NUM_BASES> pow(const MF& mf,
                                 const std::array<V, NUM_BASES>& bases,
                                 typename MF::IntegerType exponent,
                                 PTAG)
    {
        static_assert(0 < NUM_BASES && NUM_BASES <= 3, "");
        return montgomery_pow<MF>::arraypow_cond_branch_unrolled(mf, bases, exponent, PTAG());
    }
};

// Performance note: On x86-64, with TagMontyQuarterrange or TagMontyHalfrange,
// and 1 to 3 bases and compiling with clang, arraypow_cmov() always had better
// measured performance than any of the alternative functions.  For gcc,
// arraypow_cond_branch_unrolled() tended to do best.  However with gcc and 2
// bases, arraypow_masked() performed 5-10% better than
// arraypow_cond_branch_unrolled().  For simplicity, I assume this was
// anomalous.  In general we expect conditional moves (cmovs) to perform fairly
// well for the small trial sizes below, so for other compilers than clang and
// gcc we default to using arraypow_cmov() below.  In principal, the fact that
// cmov avoids conditional branches should increase the CPU's ability to exploit
// instruction level parallelism.
// Further note: On x86-64, with TagMontyFullrange (handled above by the primary
// template) and 1 to 3 bases, arraypow_cond_branch_unrolled() tended to always
// perform best.
template<class MF>
struct montgomery_array_pow_small<TagMontyHalfrange, MF> {
    using V = typename MF::MontgomeryValue;

    template <std::size_t NUM_BASES, class PTAG>
    static HURCHALLA_FORCE_INLINE
    std::array<V, NUM_BASES> pow(const MF& mf,
                                 const std::array<V, NUM_BASES>& bases,
                                 typename MF::IntegerType exponent,
                                 PTAG)
    {
        static_assert(0 < NUM_BASES && NUM_BASES <= 3, "");
#if !defined(HURCHALLA_AVOID_CSELECT) && \
    (!defined(__GNUC__) || defined(__clang__))
        return montgomery_pow<MF>::arraypow_cmov(mf, bases, exponent, PTAG());
#else
        return montgomery_pow<MF>::arraypow_cond_branch_unrolled(mf, bases, exponent, PTAG());
#endif
    }
};
template<class MF>
struct montgomery_array_pow_small<TagMontyQuarterrange, MF> {
    using V = typename MF::MontgomeryValue;

    template <std::size_t NUM_BASES, class PTAG>
    static HURCHALLA_FORCE_INLINE
    std::array<V, NUM_BASES> pow(const MF& mf,
                                 const std::array<V, NUM_BASES>& bases,
                                 typename MF::IntegerType exponent,
                                 PTAG)
    {
        static_assert(0 < NUM_BASES && NUM_BASES <= 3, "");
#if !defined(HURCHALLA_AVOID_CSELECT) && \
    (!defined(__GNUC__) || defined(__clang__))
        return montgomery_pow<MF>::arraypow_cmov(mf, bases, exponent, PTAG());
#else
        return montgomery_pow<MF>::arraypow_cond_branch_unrolled(mf, bases, exponent, PTAG());
#endif
    }
};



// Primary template.
// MF should be a MontgomeryForm type.
// The template parameter Enable should never be explicitly specified (just use
// the default).  Its purpose is to partially specialize this struct for when
// MF::IntegerType has a bit size larger than HURCHALLA_TARGET_BIT_WIDTH.
template<class MontyTag, class MF, class Enable = void>
struct montgomery_array_pow {
    using T = typename MF::IntegerType;
    using V = typename MF::MontgomeryValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");

    template <std::size_t NUM_BASES>
    static
    typename std::enable_if<(NUM_BASES <= 3), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, const std::array<V, NUM_BASES>& bases, T exponent)
    {
        static_assert(NUM_BASES > 0, "");
        using PTAG = typename std::conditional<(NUM_BASES < 2), LowlatencyTag, LowuopsTag>::type;
        return montgomery_array_pow_small<MontyTag, MF>::pow(mf, bases, exponent, PTAG());
    }

    // When we have a lot of bases we prefer to use arraypow_cond_branch()
    // instead of arraypow_cond_branch_unrolled(), to keep code size down and to
    // avoid wasting i-cache.  Arraypow_cond_branch_unrolled() force-unrolls
    // loops, which is bad for code size and i-cache, but unrolling sometimes
    // improves performance with a small NUM_BASES.
    // Note also: clang in debug mode can use up a very large amount of memory
    // for arraypow_cond_branch_unrolled() if the NUM_BASES parameter is big,
    // (enough to crash a system while testing sometimes) during compilation of
    // the Unroll class.  Since the plain arraypow_cond_branch() doesn't use
    // Unroll, it has no problem in clang when NUM_BASES starts to get big.
    //
    // The exact cutoff point between those two functions at NUM_BASES <= 5 is a
    // bit arbitrary, but probably fine.  Fyi, we would usually expect
    // 1 <= NUM_BASES <= 3, or less commonly 4 <= NUM_BASES < 8.  Otherwise
    // (much less common), we might have a "huge" NUM_BASES (perhaps >50).
    template <std::size_t NUM_BASES>
    static
    typename std::enable_if<(3 < NUM_BASES) && (NUM_BASES <= 5),
                             std::array<V, NUM_BASES>>::type
    pow(const MF& mf, const std::array<V, NUM_BASES>& bases, T exponent)
    {
        return montgomery_pow<MF>::arraypow_cond_branch_unrolled(mf, bases, exponent, LowuopsTag());
    }

    template <std::size_t NUM_BASES>
    static
    typename std::enable_if<(NUM_BASES > 5), std::array<V, NUM_BASES>>::type
    pow(const MF& mf, const std::array<V, NUM_BASES>& bases, T exponent)
    {
        return montgomery_pow<MF>::arraypow_cond_branch(mf, bases, exponent, LowuopsTag());
    }
};

// partial specialization for very large MF::IntegerTypes (e.g. __uint128_t)
template<class MontyTag, class MF>
struct montgomery_array_pow<MontyTag, MF, typename std::enable_if<
                    (ut_numeric_limits<typename MF::IntegerType>::digits
                     > HURCHALLA_TARGET_BIT_WIDTH)>::type> {
    using T = typename MF::IntegerType;
    using V = typename MF::MontgomeryValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");

    template <std::size_t NUM_BASES>
    static std::array<V, NUM_BASES>
    pow(const MF& mf, const std::array<V, NUM_BASES>& bases, T exponent)
    {
        // Note: on x86-64, 128 bit types T consistently performed best when
        // using simple arraypow_cond_branch().
        // Whenever type T is larger than the native registers, the operations
        // needed for arithmetic on T will likely expose some instruction level
        // parallelism, and so we might expect arraypow_cond_branch() to be
        // more favorable here than it would be with small T.
        using PTAG = typename std::conditional<(NUM_BASES < 2), LowlatencyTag, LowuopsTag>::type;
        return montgomery_pow<MF>::arraypow_cond_branch(mf, bases, exponent, PTAG());
    }
};


}} // end namespace

#endif
