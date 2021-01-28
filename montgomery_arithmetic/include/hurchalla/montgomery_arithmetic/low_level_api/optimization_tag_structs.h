// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

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

struct LowlatencyTag : public PrivateAnyTag {};
struct LowuopsTag : public PrivateAnyTag {};


} // end namespace


#endif
