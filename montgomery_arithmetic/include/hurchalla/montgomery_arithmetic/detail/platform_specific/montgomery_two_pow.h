// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_TWO_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_TWO_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/impl_montgomery_two_pow.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyTags.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstddef>
#include <array>
#include <type_traits>

namespace hurchalla { namespace detail {


class Tag_montgomery_two_pow_gcc {};
class Tag_montgomery_two_pow_clang {};

// Essentially Tag_montgomery_two_pow_big is for __uint128_t (or bigger?)
class Tag_montgomery_two_pow_big {};
// Essentially Tag_montgomery_two_pow_small is for uint64_t or smaller
class Tag_montgomery_two_pow_small {};



// Primary template for tagging different perf tunings
template <class MontyTag, class CompilerTag, class PowSizeTag>
struct tagged_montgomery_two_pow {};




// These specializations call impl_montgomery_two_pow, using the template
// arguments that were found to perform best in benchmarks for a corresponding
// configuration of compiler/uint_type/Mont_type


// Fyi, the meaning of the impl_montgomery_two_pow::call template parameters are
// scalar call
//  template <class MF, typename U,
//            bool USE_SLIDING_WINDOW_OPTIMIZATION,
//            size_t TABLE_BITS, size_t CODE_SECTION,
//            bool USE_SQUARING_VALUE_OPTIMIZATION>
//
// array call
//  template <class MF, typename U,
//            size_t ARRAY_SIZE, size_t TABLE_BITS, size_t CODE_SECTION,
//            bool USE_SQUARING_VALUE_OPTIMIZATION>


// -- the following best performance tunings were measured with mac M2 --


// Partial Specialization: clang and big uint pow.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // false 0 34 true probably best.
    // false 0 33 true best for asm, but with noasm it is ~10% slower than
    // 34.  Whereas with asm, 34 is ~1% slower than 33.  Although I strongly
    // favor asm results, best compromise is using 34.
    return impl_montgomery_two_pow::call<MF, U, false, 0, 34, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 31 best.  30 trails ~0.15% with asm, and trails ~0.6% with no asm
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: clang and big uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // false 0 34 false clear winner
    return impl_montgomery_two_pow::call<MF, U, false, 0, 34, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // slightly prefer 0 31 false due to tiny bit better on asm.  Very close to a
    // toss-up with 0 30 false.
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: clang and big uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // False 0 34 false and 22 are equally good (34 wins on asm by ~0.15% and
    // loses by ~0.7% on noasm)
    return impl_montgomery_two_pow::call<MF, U, false, 0, 34, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 0, 31, false best.  It's ~2.5% faster than 30 for asm (and 1% faster
    // than 29), though ~1.5% slower than 30 with noasm.
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};


// Partial specialization: gcc and big uint pow.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // prefer false 0 22 true.  false 0 34 true a very good second choice.
    // false 0 22 true best for asm, and ~1% slower than 0 34 true on
    // noasm (~2.5% slower than noasm best).  false 0 34 true a good compromise
    // to help noasm- for asm, it is ~0.5% slower than asm best.
    return impl_montgomery_two_pow::call<MF, U, false, 0, 22, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 0 31 false prefered.  0 30 false trails by ~0.4% with asm, and leads by
    // ~0.25% with noasm.  Note that for noasm, 0 30 true leads by ~4.5% over
    // 0 30 false; however using true loses ~6% perf with asm, so it would be a
    // bad compromise.
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: gcc and big uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // prefer false 0 22 false.  34 is also good - it is ~0.4% slower with asm
    // and ~0.6% slower with noasm
    return impl_montgomery_two_pow::call<MF, U, false, 0, 22, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // Use 0 31 false:
    // We slightly prefer 31 for asm (very close to toss-up with 30), and 31 is
    // ~2% faster than 30 with noasm.
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: gcc and big uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // False 0 22 false and 34 are equally good.  34 is ~0.15% slower than 22
    // on asm but 0.45% faster than 22 with noasm.
    return impl_montgomery_two_pow::call<MF, U, false, 0, 22, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 0 31 false wins.  It's a tiny bit better on asm than 30 (almost a toss-
    // up) and very roughly over 1% faster than 30 on noasm.
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};



// Partial Specialization: clang and small uint pow.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // false, 0, 21, true  best, due to asm.  Second place is
    // false, 0, 23, true - it's a bit of a compromise of asm/noasm - it's 0.45%
    // slower than 21 at asm, and 0.45% faster at noasm.
    //
    // For consistency we'll choose 23 since it's used elsewhere and 21 is not.
    return impl_montgomery_two_pow::call<MF, U, false, 0, 23, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 0 31 false best, both asm, noasm
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: clang and small uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // 25 or 37 almost equal.  prefer a tiny amount 25 (it wins by ~0.1% at asm
    // and loses ~0.1% at noasm).  35 would be next best choice
    return impl_montgomery_two_pow::call<MF, U, false, 0, 25, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 31 preferred, but 29 wins at biggest array sizes by ~0.25% - 0.5%.  So 29
    // is fine too (it loses by larger margins at smaller sizes).
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 29, false>(mf, n);
  }
};
// Full Specialization: clang and small uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // 23 best, 35 next best and is 0.2-0.4% slower
    return impl_montgomery_two_pow::call<MF, U, false, 0, 23, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 31 wins overall, 29 very close and arguably wins with noasm.
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};


// Partial Specialization: gcc and small uint pow and MontgomeryFull.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // false, 0, 23, true wins.  this prioritizes asm, and is best compromise
    // with noasm too.
    return impl_montgomery_two_pow::call<MF, U, false, 0, 23, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 0 29 false and 31 equally good.  for absolute perf 31 wins by 0.1 - 0.5%,
    // at cost of biggest array sizes.  29 wins at smaller sizes by a larger
    // margin (around 0.5% - 1%, and up to 2%).
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: gcc and small uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // 25 best, though 37 behind only 0.1 - 0.25%
    return impl_montgomery_two_pow::call<MF, U, false, 0, 25, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 0, 28, false wins.  30 fine too though a bit slower
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 28, false>(mf, n);
  }
};
// Full Specialization: gcc and small uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    // 37 wins, 25 second place - 0.2% slower.
    return impl_montgomery_two_pow::call<MF, U, false, 0, 37, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    // 29 seems best (it wins at biggest array sizes), but it's not completely
    // clear what's best overall.  28 seems best overall if you favor small
    // array sizes.
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 29, false>(mf, n);
  }
};




struct montgomery_two_pow {

  // Calculate pow(2, n), modulo the modulus of mf, and return the result in
  // montgomeryform representation.
  template <class MF, typename T>
  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, T nt)
  {
    HPBC_CLOCKWORK_PRECONDITION(nt >= 0);
    using U = typename extensible_make_unsigned<T>::type;
    U n = static_cast<U>(nt);

    using MontyTag = typename MF::MontType::MontyTag;
    using RU = typename MontgomeryFormExtensions<MF, LowlatencyTag>::RU;
    constexpr bool isBigPow = ut_numeric_limits<RU>::digits >
                              HURCHALLA_TARGET_BIT_WIDTH;
    using SizeTag = typename std::conditional<isBigPow,
                Tag_montgomery_two_pow_big, Tag_montgomery_two_pow_small>::type;
#if (!defined(__GNUC__) || defined(__clang__))
    using Compiler = Tag_montgomery_two_pow_clang;   // clang tuning
#else
    using Compiler = Tag_montgomery_two_pow_gcc;     // gcc tuning
#endif
    
    return tagged_montgomery_two_pow<MontyTag, Compiler, SizeTag>::call(mf, n);
  }


  // Helper function - delegated Array version of montgomery two pow
  template <class MF, typename U, std::size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  helper(const std::array<MF,ARRAY_SIZE>& mf, const std::array<U,ARRAY_SIZE>& n)
  {
    static_assert(hurchalla::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hurchalla::ut_numeric_limits<U>::is_signed, "");

    using MontyTag = typename MF::MontType::MontyTag;
    using RU = typename MontgomeryFormExtensions<MF, LowlatencyTag>::RU;
    constexpr bool isBigPow = ut_numeric_limits<RU>::digits >
                              HURCHALLA_TARGET_BIT_WIDTH;
    using SizeTag = typename std::conditional<isBigPow,
                Tag_montgomery_two_pow_big, Tag_montgomery_two_pow_small>::type;
#if (!defined(__GNUC__) || defined(__clang__))
    using Compiler = Tag_montgomery_two_pow_clang;   // clang tuning
#else
    using Compiler = Tag_montgomery_two_pow_gcc;     // gcc tuning
#endif

    return tagged_montgomery_two_pow<MontyTag, Compiler, SizeTag>::call(mf, n);
  }


  // Array version of montgomery two pow, for unsigned T
  template <class MF, typename T, std::size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE static
  typename std::enable_if<!(hurchalla::ut_numeric_limits<T>::is_signed),
                          std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
                          >::type
  call(const std::array<MF,ARRAY_SIZE>& mf, const std::array<T,ARRAY_SIZE>& nt)
  {
      return helper(mf, nt);
  }

  // Array version of montgomery two pow, for signed T
  template <class MF, typename T, std::size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE static
  typename std::enable_if<(hurchalla::ut_numeric_limits<T>::is_signed),
                          std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
                          >::type
  call(const std::array<MF,ARRAY_SIZE>& mf, const std::array<T,ARRAY_SIZE>& nt)
  {
      using U = typename extensible_make_unsigned<T>::type;
      std::array<U, ARRAY_SIZE> n;
      HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<ARRAY_SIZE; ++i) {
        HPBC_CLOCKWORK_PRECONDITION(nt[i] >= 0);
        n[i] = static_cast<U>(nt[i]);
      }
      return helper(mf, n);
  }
};


}} // end namespace

#endif
