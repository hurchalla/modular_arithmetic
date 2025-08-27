// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_PROGRAMMING_BY_CONTRACT_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_PROGRAMMING_BY_CONTRACT_H_INCLUDED


#ifdef _MSC_VER
#  define HPBC_CLOCKWORK_DO_NOTHING(...) (false ? (void)(__VA_ARGS__) : (void)0)
#else
#  define HPBC_CLOCKWORK_DO_NOTHING(...) ((void)(true || (__VA_ARGS__)))
#endif



#if defined(HURCHALLA_CLOCKWORK_ENABLE_ASSERTS) || \
    defined(CLOCKWORK_CHECK_API_PRECONDITIONS)

#  include <cstdio>
#  include <cstdlib>
      /* We can (probably) detect if exceptions are enabled by checking the
         gcc/clang macro __EXCEPTIONS, msvc's macro _CPPUNWIND, and the official
         but not always supported C++98 macro __cpp_exceptions. */
#  if defined(__cplusplus) && (defined(__EXCEPTIONS) || defined(_CPPUNWIND) \
                        || (defined(__cpp_exceptions) && __cpp_exceptions != 0))
      // with exceptions: treat an exception as a failure during assert
#     define HPBC_CLOCKWORK_BASIC_ASSERT(...) \
               do { \
                  bool assertPassed = false; \
                  try { if (__VA_ARGS__) assertPassed = true; } \
                  catch (...) {} \
                  if (!assertPassed) { \
                     fprintf(stderr, "Assert failed (%s): file %s, line %d\n", \
                              #__VA_ARGS__, __FILE__, __LINE__); \
                     std::abort(); \
                  } \
               } while(0)
#  else   /* without exceptions */
#     define HPBC_CLOCKWORK_BASIC_ASSERT(...) \
               do { \
                  if (__VA_ARGS__) {} \
                  else { \
                     fprintf(stderr, "Assert failed (%s): file %s, line %d\n", \
                              #__VA_ARGS__, __FILE__, __LINE__); \
                     std::abort(); \
                  } \
               } while(0)
#  endif

#  define HPBC_CLOCKWORK_API_PRECONDITION(...) do { \
                            HPBC_CLOCKWORK_BASIC_ASSERT(__VA_ARGS__); } while(0)

#else
#  define HPBC_CLOCKWORK_API_PRECONDITION(...) \
                            HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)

#endif



#if !defined(HURCHALLA_CLOCKWORK_ENABLE_ASSERTS)

#  define HPBC_CLOCKWORK_PRECONDITION(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_PRECONDITION2(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_PRECONDITION3(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_POSTCONDITION(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_POSTCONDITION2(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_POSTCONDITION3(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_INVARIANT(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_INVARIANT2(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_INVARIANT3(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_ASSERT(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_ASSERT2(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_ASSERT3(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)

#  if defined(__cplusplus)
#     define HPBC_CLOCKWORK_FALSE_VALUE (false)
#  else
#     define HPBC_CLOCKWORK_FALSE_VALUE (0)
#  endif
#  define HPBC_CLOCKWORK_PRECONDITION_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_PRECONDITION3_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_POSTCONDITION_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_POSTCONDITION3_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_INVARIANT_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_INVARIANT2_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_INVARIANT3_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_ASSERT_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_ASSERT2_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE
#  define HPBC_CLOCKWORK_ASSERT3_MACRO_IS_ACTIVE HPBC_CLOCKWORK_FALSE_VALUE

#else

#  define HPBC_CLOCKWORK_DEFAULT_ASSERT_LEVEL       3

#  if defined(HURCHALLA_CLOCKWORK_ASSERT_LEVEL) && (1 - HURCHALLA_CLOCKWORK_ASSERT_LEVEL - 1 == 2)
#    error "HURCHALLA_CLOCKWORK_ASSERT_LEVEL is defined but with no value, or negative value"
#  endif
#  if !defined(HURCHALLA_CLOCKWORK_ASSERT_LEVEL)
#    define HURCHALLA_CLOCKWORK_ASSERT_LEVEL  HPBC_CLOCKWORK_DEFAULT_ASSERT_LEVEL
#  endif
#  if ((HURCHALLA_CLOCKWORK_ASSERT_LEVEL < 0) || (HURCHALLA_CLOCKWORK_ASSERT_LEVEL > 3))
#    error "Invalid assert level for HURCHALLA_CLOCKWORK_ASSERT_LEVEL"
#  endif

#  define HPBC_CLOCKWORK_LEVEL_ASSERT(LEVEL, ...) do { \
                       if (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= LEVEL) { \
                            HPBC_CLOCKWORK_BASIC_ASSERT(__VA_ARGS__); } } while(0)

#  define HPBC_CLOCKWORK_PRECONDITION(...) HPBC_CLOCKWORK_LEVEL_ASSERT(1, __VA_ARGS__)
#  define HPBC_CLOCKWORK_PRECONDITION2(...) HPBC_CLOCKWORK_LEVEL_ASSERT(2, __VA_ARGS__)
#  define HPBC_CLOCKWORK_PRECONDITION3(...) HPBC_CLOCKWORK_LEVEL_ASSERT(3, __VA_ARGS__)
#  define HPBC_CLOCKWORK_POSTCONDITION(...) HPBC_CLOCKWORK_LEVEL_ASSERT(1, __VA_ARGS__)
#  define HPBC_CLOCKWORK_POSTCONDITION2(...) HPBC_CLOCKWORK_LEVEL_ASSERT(2, __VA_ARGS__)
#  define HPBC_CLOCKWORK_POSTCONDITION3(...) HPBC_CLOCKWORK_LEVEL_ASSERT(3, __VA_ARGS__)
#  define HPBC_CLOCKWORK_INVARIANT(...) HPBC_CLOCKWORK_LEVEL_ASSERT(1, __VA_ARGS__)
#  define HPBC_CLOCKWORK_INVARIANT2(...) HPBC_CLOCKWORK_LEVEL_ASSERT(2, __VA_ARGS__)
#  define HPBC_CLOCKWORK_INVARIANT3(...) HPBC_CLOCKWORK_LEVEL_ASSERT(3, __VA_ARGS__)
#  define HPBC_CLOCKWORK_ASSERT(...) HPBC_CLOCKWORK_LEVEL_ASSERT(1, __VA_ARGS__)
#  define HPBC_CLOCKWORK_ASSERT2(...) HPBC_CLOCKWORK_LEVEL_ASSERT(2, __VA_ARGS__)
#  define HPBC_CLOCKWORK_ASSERT3(...) HPBC_CLOCKWORK_LEVEL_ASSERT(3, __VA_ARGS__)

#  define HPBC_CLOCKWORK_PRECONDITION_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 1)
#  define HPBC_CLOCKWORK_PRECONDITION2_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 2)
#  define HPBC_CLOCKWORK_PRECONDITION3_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 3)
#  define HPBC_CLOCKWORK_POSTCONDITION_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 1)
#  define HPBC_CLOCKWORK_POSTCONDITION2_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 2)
#  define HPBC_CLOCKWORK_POSTCONDITION3_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 3)
#  define HPBC_CLOCKWORK_INVARIANT_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 1)
#  define HPBC_CLOCKWORK_INVARIANT2_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 2)
#  define HPBC_CLOCKWORK_INVARIANT3_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 3)
#  define HPBC_CLOCKWORK_ASSERT_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 1)
#  define HPBC_CLOCKWORK_ASSERT2_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 2)
#  define HPBC_CLOCKWORK_ASSERT3_MACRO_IS_ACTIVE (HURCHALLA_CLOCKWORK_ASSERT_LEVEL >= 3)

#endif



#if defined(HURCHALLA_CLOCKWORK_ENABLE_ASSERTS)
   // this section was adapted from the ideas in
   // https://akrzemi1.wordpress.com/2017/05/18/asserts-in-constexpr-functions/
   // https://gist.github.com/oliora/928424f7675d58fadf49c70fdba70d2f
#  include "hurchalla/util/compiler_macros.h"
#  include <utility>
#  if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#  endif
   template <class L>
   void hurchalla_hpbc_clockwork_forward_lambda(L&& lambda) noexcept 
   {
      std::forward<L>(lambda)();
   }
#  if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic pop
#  endif
#  ifdef _MSC_VER
#    define HPBC_CLOCKWORK_CONSTEXPR_ASSERT(...) ((void)(HURCHALLA_LIKELY(__VA_ARGS__) ? \
                             (void)0 : hurchalla_hpbc_clockwork_forward_lambda( \
                             [](){ HPBC_CLOCKWORK_ASSERT(!#__VA_ARGS__);}), (void)0))
#  else
#    define HPBC_CLOCKWORK_CONSTEXPR_ASSERT(...) ((void)(HURCHALLA_LIKELY(__VA_ARGS__) ? \
                             (void)0 : hurchalla_hpbc_clockwork_forward_lambda( \
                             [](){ HPBC_CLOCKWORK_ASSERT(#__VA_ARGS__ == nullptr);}), (void)0))
#  endif
#  define HPBC_CLOCKWORK_CONSTEXPR_PRECONDITION(...) HPBC_CLOCKWORK_CONSTEXPR_ASSERT(__VA_ARGS__)
#  define HPBC_CLOCKWORK_CONSTEXPR_POSTCONDITION(...) HPBC_CLOCKWORK_CONSTEXPR_ASSERT(__VA_ARGS__)
#  define HPBC_CLOCKWORK_CONSTEXPR_INVARIANT(...) HPBC_CLOCKWORK_CONSTEXPR_ASSERT(__VA_ARGS__)

#else

#  define HPBC_CLOCKWORK_CONSTEXPR_ASSERT(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_CONSTEXPR_PRECONDITION(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_CONSTEXPR_POSTCONDITION(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)
#  define HPBC_CLOCKWORK_CONSTEXPR_INVARIANT(...) HPBC_CLOCKWORK_DO_NOTHING(__VA_ARGS__)

#endif


#endif
