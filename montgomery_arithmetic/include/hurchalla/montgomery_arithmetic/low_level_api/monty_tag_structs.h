// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_TAG_STRUCTS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_TAG_STRUCTS_H_INCLUDED


namespace hurchalla {


struct FullrangeTag {};
struct QuarterrangeTag {};

// The name "Fullrange" signifies that there are essentially no preconditions on
// the value of the modulus.  Although montgomery multiplication always requires
// that the modulus is odd, this function will work for any (odd) modulus that
// is representable by its type T.
//
// The name "Quarterrange" signifies that the modulus must be less than R/4,
// where  R = 2^(ut_numeric_limits<T>::digits).  For example, if T is uint64_t
// then R = 2^64 and R/4 == 2^62, and thus it would require  modulus < 2^62.
//
// Quarterrange functions require/allow an unusual input range:
// for an input x, they allow  0 <= x < 2*n, where n is the modulus.  Similarly,
// the return value range will be  0 <= returnValue < 2*n.  Obviously neither
// inputs nor outputs necessarily belong to the minimal residue class modulo n -
// i.e. they might not be fully reduced, modulo n.  Note that the algorithm for
// montgomery REDC requires that  u = x*y < n*R;  this will always be satisfied
// for any multiplication x*y of Quarterrange montgomery values.  To see why,
// keep in mind that Quarterrange requires n < R/4  and that all inputs are less
// than 2*n.  Thus the multiplication
// u = x*y < (2*n)*(2*n) == (4*n)*n < (4*n)*(R/4) == n*R,  which means u < n*R,
// as required.  For more details on Quarterrange, see also section 5 from
// "Montgomery's Multiplication Technique: How to Make It Smaller and Faster"
// https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps


} // end namespace


#endif
