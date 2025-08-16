// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_TAGS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_TAGS_H_INCLUDED



namespace hurchalla { namespace detail {


struct TagMontyQuarterrange final {};  // IDs MontyQuarterRange independent of T
struct TagMontyHalfrange final {};
struct TagMontyFullrange final {};
struct TagMontyWrappedmath final {};
struct TagMontyFullrangeMasked final {};


}} // end namespace

#endif
