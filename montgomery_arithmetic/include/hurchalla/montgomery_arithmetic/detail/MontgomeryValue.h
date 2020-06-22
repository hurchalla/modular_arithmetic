// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_VALUE_H_INCLUDED


#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Weffc++"
#endif


namespace hurchalla { namespace montgomery_arithmetic {


// A simple wrapper for T, used to designate that a value is in Montgomery Form
template<typename T> struct MontgomeryValue final {
public:
    MontgomeryValue() {} // This constructor purposely does not initialize
                                   // 'value' - the contents are undefined
                                   // until the object is assigned to.
    explicit MontgomeryValue(T val) : value(val) {}
    T get() const { return value; }
private:
    T value;
};

template<typename T>
bool operator==(const MontgomeryValue<T>& lhs, const MontgomeryValue<T>& rhs) {
    return lhs.get() == rhs.get();
}
template<typename T>
bool operator!=(const MontgomeryValue<T>& lhs, const MontgomeryValue<T>& rhs) {
    return !(lhs == rhs);
}


}} // end namespace


#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif


#endif
