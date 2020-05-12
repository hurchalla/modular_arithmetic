
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_VALUE_H_INCLUDED

namespace hurchalla { namespace montgomery_arithmetic {


// A simple wrapper for T, used to designate that a value is in Montgomery Form
template<typename T> struct MontgomeryValue final {
public:
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

#endif
