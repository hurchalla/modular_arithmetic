// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_OPTIMIZATION_TAG_STRUCTS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_OPTIMIZATION_TAG_STRUCTS_H_INCLUDED


namespace hurchalla { namespace montgomery_arithmetic {


struct PrivateInplaceTag {};
struct PrivateOutofplaceTag {};

struct InplaceLowlatencyTag : public PrivateInplaceTag {};
struct OutofplaceLowlatencyTag : public PrivateOutofplaceTag {};
struct InplaceLowuopsTag : public PrivateInplaceTag {};
struct OutofplaceLowuopsTag : public PrivateOutofplaceTag {};


struct FullrangeTag {};
struct HalfrangeTag {};
struct QuarterrangeTag {};


}} // end namespace


#endif
