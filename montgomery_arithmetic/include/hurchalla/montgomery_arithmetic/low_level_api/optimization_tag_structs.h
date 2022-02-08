// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_OPTIMIZATION_TAG_STRUCTS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_OPTIMIZATION_TAG_STRUCTS_H_INCLUDED


namespace hurchalla {


// private optimization tag intended only for use by the implementation
struct PrivateAnyTag {};


// Public optimization tags:
// ------------------------
// LowlatencyTag potentially offers optimizations targeted toward lowering the
// latency of functions.
// LowuopsTag potentially offers optimizations targeted toward reducing the
// number of instructions generated/executed by functions.

struct LowlatencyTag final : public PrivateAnyTag {};
struct LowuopsTag final : public PrivateAnyTag {};


} // end namespace


#endif
