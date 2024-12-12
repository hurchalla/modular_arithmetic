// Copyright (c) 2024 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>
#include <cstddef>

namespace hurchalla { namespace detail {


// The primary template below handles when InlineAll == true, and annotates
// all class functions with a force inline attribute.
// The template specialization handles InlineAll == false, and does not annotate
// any of the class functionw with a force inline attibute.
//
// Note: this is a rare case where ugly #define / #undef / #include hacking
// seems to be the best way to make the code clear and maintainable.  Placing or
// not placing an attribute on a function doesn't appear to be something we can
// directly do with a template parameter.  So we work-around it by creating two
// exact class duplicates (not counting the attribute, which is #defined or
// #undef'd) by using #include, and we use a class specialization (as below)
// to determine whether or not the class's functions get the attribute defined
// and placed, or not.


#define HURCHALLA_IMF_MAYBE_FORCE_INLINE HURCHALLA_FORCE_INLINE
//
// Primary template, instantiated for InlineAll == true.
//
// All functions in this instantiation get a force inline annotation.
template <class T, bool InlineAll, class MontyType>
class ImplMontgomeryForm final {
#include "hurchalla/montgomery_arithmetic/detail/ImplMontgomeryForm.contents"
};
#undef HURCHALLA_IMF_MAYBE_FORCE_INLINE


#define HURCHALLA_IMF_MAYBE_FORCE_INLINE
//
// Specialization, instantiated for InlineAll == false.
//
// No functions will get a force inline annotation, because
// HURCHALLA_IMF_MAYBE_FORCE_INLINE is blank.
template <class T, class MontyType>
class ImplMontgomeryForm<T, false, MontyType> final {
#include "hurchalla/montgomery_arithmetic/detail/ImplMontgomeryForm.contents"
};
#undef HURCHALLA_IMF_MAYBE_FORCE_INLINE


}} // end namespace

#endif
