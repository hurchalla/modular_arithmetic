
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_VALUE_H__INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_VALUE_H__INCLUDED

namespace hurchalla { namespace montgomery_arithmetic {


// A simple wrapper for T, used to designate that a value is in Montgomery Form
template<typename T> struct MontgomeryValue final {
public:
    MontgomeryValue() {}; // This constructor purposely does not initialize
                                   // 'value' - the contents are undefined
                                   // until the object is assigned to.
    explicit MontgomeryValue(T val) : value(val) {}
    T get() { return value; }
private:
    T value;
};

template<typename T>
bool operator==(MontgomeryValue<T>& lhs, MontgomeryValue<T>& rhs) {
    return lhs.get() == rhs.get();
}
template<typename T>
bool operator!=(MontgomeryValue<T>& lhs, MontgomeryValue<T>& rhs) {
    return !(lhs == rhs);
}


}} // end namespace

#endif
