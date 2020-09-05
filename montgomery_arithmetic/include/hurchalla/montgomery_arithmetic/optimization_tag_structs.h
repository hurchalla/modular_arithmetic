// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_OPTIMIZATION_TAG_STRUCTS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_OPTIMIZATION_TAG_STRUCTS_H_INCLUDED


namespace hurchalla { namespace montgomery_arithmetic {


// private optimization tags intended only for use by the implementation
struct PrivateAnyTag {};
struct PrivateInplaceTag : public PrivateAnyTag {};
struct PrivateOutofplaceTag : public PrivateAnyTag {};

// public optimization tags
struct InplaceLowlatencyTag : public PrivateInplaceTag {};
struct OutofplaceLowlatencyTag : public PrivateOutofplaceTag {};
struct InplaceLowuopsTag : public PrivateInplaceTag {};
struct OutofplaceLowuopsTag : public PrivateOutofplaceTag {};


}} // end namespace


#endif
