INT_RUNTIME = r'''#include <boost/multiprecision/cpp_int.hpp>
#include <memory>
#include <vector>
#include <limits>
#include <cmath>
#include <cctype>
#include <compare>
#include <string>
#include <iostream>
#include <type_traits>
#include <stdexcept>
#include <functional>
#include <utility>

class Int {
public:
    using Big = boost::multiprecision::cpp_int;
private:
    long long small_ = 0;
    std::unique_ptr<Big> big_;
    explicit Int(Big value) {
        if (value >= std::numeric_limits<long long>::min() && value <= std::numeric_limits<long long>::max()) small_ = value.convert_to<long long>();
        else big_ = std::make_unique<Big>(std::move(value));
    }
    static Int trunc_div(const Int& a, const Int& b) {
        if (b == 0) throw std::domain_error("integer division by zero");
        if (!a.big_ && !b.big_ && !(a.small_ == std::numeric_limits<long long>::min() && b.small_ == -1)) return Int(a.small_ / b.small_);
        return Int(Big(a.big_value() / b.big_value()));
    }
public:
    Int() = default;
    template<class T> requires std::is_integral_v<T>
    Int(T value) {
        if constexpr (std::is_unsigned_v<T>) {
            if (value > static_cast<unsigned long long>(std::numeric_limits<long long>::max())) { big_ = std::make_unique<Big>(value); return; }
        }
        if constexpr (sizeof(T) > sizeof(long long) && std::is_signed_v<T>) {
            if (value < std::numeric_limits<long long>::min() || value > std::numeric_limits<long long>::max()) { big_ = std::make_unique<Big>(value); return; }
        }
        small_ = static_cast<long long>(value);
    }
    Int(const Int& value) : small_(value.small_), big_(value.big_ ? std::make_unique<Big>(*value.big_) : nullptr) {}
    Int(Int&&) noexcept = default;
    Int& operator=(const Int& value) { if (this != &value) { small_ = value.small_; big_ = value.big_ ? std::make_unique<Big>(*value.big_) : nullptr; } return *this; }
    Int& operator=(Int&&) noexcept = default;
    template<class T> requires std::is_floating_point_v<T>
    explicit Int(T value) {
        if (!std::isfinite(value)) throw std::domain_error("cannot convert nonfinite float to integer");
        *this = Int(Big(value));
    }
    static Int from_big(Big value) { return Int(std::move(value)); }
    Big big_value() const { return big_ ? *big_ : Big(small_); }
    bool is_small() const { return !big_; }
    long long small_value() const { if (big_) throw std::overflow_error("integer does not fit in int64"); return small_; }
    explicit operator long long() const { return small_value(); }
    template<class T> requires (std::is_integral_v<T> && !std::is_same_v<T,bool>)
    operator T() const {
        if (!big_) {
            if (std::in_range<T>(small_)) return static_cast<T>(small_);
            throw std::overflow_error("integer conversion out of range");
        }
        if (*big_ < std::numeric_limits<T>::min() || *big_ > std::numeric_limits<T>::max()) throw std::overflow_error("integer conversion out of range");
        return big_->convert_to<T>();
    }
    template<class T> requires std::is_floating_point_v<T>
    explicit operator T() const { return big_ ? big_->convert_to<T>() : static_cast<T>(small_); }
    explicit operator bool() const { return big_ ? *big_ != 0 : small_ != 0; }
    static Int parse(const std::string& text, int base = 10) {
        size_t pos = 0;
        while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
        bool negative = false;
        if (pos < text.size() && (text[pos] == '-' || text[pos] == '+')) negative = text[pos++] == '-';
        if (base == 0) {
            base = 10;
            if (pos + 1 < text.size() && text[pos] == '0') {
                char marker = text[pos + 1];
                if (marker == 'x' || marker == 'X') base = 16;
                else if (marker == 'b' || marker == 'B') base = 2;
                else if (marker == 'o' || marker == 'O') base = 8;
            }
        }
        if (base < 2 || base > 36) throw std::invalid_argument("invalid integer base");
        if (pos + 1 < text.size() && text[pos] == '0' && ((base == 16 && (text[pos+1] == 'x' || text[pos+1] == 'X')) || (base == 2 && (text[pos+1] == 'b' || text[pos+1] == 'B')) || (base == 8 && (text[pos+1] == 'o' || text[pos+1] == 'O')))) pos += 2;
        Int result = 0;
        bool digit_seen = false, underscore = false;
        for (; pos < text.size(); ++pos) {
            unsigned char c = text[pos];
            if (std::isspace(c)) break;
            if (c == '_') { if (!digit_seen || underscore) throw std::invalid_argument("invalid integer"); underscore = true; continue; }
            int digit = c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'z' ? c - 'a' + 10 : c >= 'A' && c <= 'Z' ? c - 'A' + 10 : -1;
            if (digit < 0 || digit >= base) throw std::invalid_argument("invalid integer");
            result = result * base + digit; digit_seen = true; underscore = false;
        }
        if (!digit_seen || underscore) throw std::invalid_argument("invalid integer");
        while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
        if (pos != text.size()) throw std::invalid_argument("invalid integer");
        return negative ? -result : result;
    }
    std::string str() const { return big_ ? big_->str() : std::to_string(small_); }
    friend std::ostream& operator<<(std::ostream& out, const Int& value) { return out << value.str(); }
    friend std::istream& operator>>(std::istream& in, Int& value) { std::string token; if (in >> token) value = parse(token); return in; }
    friend bool operator==(const Int& a, const Int& b) { return !a.big_ && !b.big_ ? a.small_ == b.small_ : a.big_value() == b.big_value(); }
    friend std::strong_ordering operator<=>(const Int& a, const Int& b) {
        if (!a.big_ && !b.big_) return a.small_ <=> b.small_;
        Big left = a.big_value(), right = b.big_value(); return left < right ? std::strong_ordering::less : left > right ? std::strong_ordering::greater : std::strong_ordering::equal;
    }
    friend Int operator+(const Int& a, const Int& b) { long long result; if (!a.big_ && !b.big_ && !__builtin_add_overflow(a.small_, b.small_, &result)) return Int(result); return Int(Big(a.big_value() + b.big_value())); }
    friend Int operator-(const Int& a, const Int& b) { long long result; if (!a.big_ && !b.big_ && !__builtin_sub_overflow(a.small_, b.small_, &result)) return Int(result); return Int(Big(a.big_value() - b.big_value())); }
    friend Int operator*(const Int& a, const Int& b) { long long result; if (!a.big_ && !b.big_ && !__builtin_mul_overflow(a.small_, b.small_, &result)) return Int(result); return Int(Big(a.big_value() * b.big_value())); }
    friend Int operator/(const Int& a, const Int& b) {
        if (b == 0) throw std::domain_error("integer division by zero");
        if (!a.big_ && !b.big_ && !(a.small_ == std::numeric_limits<long long>::min() && b.small_ == -1)) {
            long long q=a.small_/b.small_, r=a.small_%b.small_;
            if (r!=0 && ((r<0)!=(b.small_<0))) --q;
            return Int(q);
        }
        Int q = trunc_div(a,b); Int r = a - q*b; if (r != 0 && ((r < 0) != (b < 0))) --q; return q; }
    friend Int operator%(const Int& a, const Int& b) { if (b == 0) throw std::domain_error("integer modulo by zero"); if (!a.big_ && !b.big_) { long long r = b.small_ == -1 ? 0 : a.small_ % b.small_; if (r != 0 && ((r < 0) != (b.small_ < 0))) r += b.small_; return Int(r); } Big r = a.big_value() % b.big_value(); if (r != 0 && ((r < 0) != (b < 0))) r += b.big_value(); return Int(std::move(r)); }
    Int operator-() const { if (!big_ && small_ != std::numeric_limits<long long>::min()) return Int(-small_); return Int(Big(-big_value())); }
    Int operator+() const { return *this; }
    Int operator~() const { if (!big_) return Int(~small_); return -(*this) - 1; }
    friend Int operator&(const Int& a, const Int& b) { if (!a.big_ && !b.big_) return Int(a.small_ & b.small_); return Int(Big(a.big_value() & b.big_value())); }
    friend Int operator|(const Int& a, const Int& b) { if (!a.big_ && !b.big_) return Int(a.small_ | b.small_); return Int(Big(a.big_value() | b.big_value())); }
    friend Int operator^(const Int& a, const Int& b) { if (!a.big_ && !b.big_) return Int(a.small_ ^ b.small_); return Int(Big(a.big_value() ^ b.big_value())); }
    friend Int operator<<(const Int& a, const Int& b) {
        if (b < 0) throw std::domain_error("negative shift count");
        if (a == 0) return 0;
        auto shift = static_cast<size_t>(b);
        if (!a.big_ && shift < 63) { long long result; if (!__builtin_mul_overflow(a.small_, (1LL << shift), &result)) return Int(result); }
        return Int(Big(a.big_value() << shift));
    }
    friend Int operator>>(const Int& a, const Int& b) {
        if (b < 0) throw std::domain_error("negative shift count");
        if (b > a.bit_length()) return a < 0 ? Int(-1) : Int(0);
        auto shift = static_cast<size_t>(b);
        if (!a.big_) return shift >= 64 ? Int(a.small_ < 0 ? -1 : 0) : Int(a.small_ >> shift);
        if (a < 0) return -Int(Big(((-a).big_value() + (Big(1) << shift) - 1) >> shift));
        return Int(Big(a.big_value() >> shift));
    }
    Int& operator++() { return *this += 1; } Int operator++(int) { Int old = *this; ++*this; return old; }
    Int& operator--() { return *this -= 1; } Int operator--(int) { Int old = *this; --*this; return old; }
    Int& operator+=(const Int& b) { return *this = *this + b; } Int& operator-=(const Int& b) { return *this = *this - b; }
    Int& operator*=(const Int& b) { return *this = *this * b; } Int& operator/=(const Int& b) { return *this = *this / b; }
    Int& operator%=(const Int& b) { return *this = *this % b; } Int& operator&=(const Int& b) { return *this = *this & b; }
    Int& operator|=(const Int& b) { return *this = *this | b; } Int& operator^=(const Int& b) { return *this = *this ^ b; }
    Int& operator<<=(const Int& b) { return *this = *this << b; } Int& operator>>=(const Int& b) { return *this = *this >> b; }
    Int isqrt() const { if (*this < 0) throw std::domain_error("isqrt() argument must be nonnegative"); return Int(Big(boost::multiprecision::sqrt(big_value()))); }
    Int bit_length() const { if (!big_) { unsigned long long magnitude = small_ < 0 ? 0ULL-static_cast<unsigned long long>(small_) : static_cast<unsigned long long>(small_); return magnitude == 0 ? 0 : 64 - __builtin_clzll(magnitude); } Big value = big_value(); if (value < 0) value = -value; return value == 0 ? Int(0) : Int(boost::multiprecision::msb(value) + 1); }
    Int bit_count() const { if (!big_) { unsigned long long magnitude = small_ < 0 ? 0ULL-static_cast<unsigned long long>(small_) : static_cast<unsigned long long>(small_); return __builtin_popcountll(magnitude); } Big value = big_value(); if (value < 0) value = -value; std::vector<unsigned long long> limbs; boost::multiprecision::export_bits(value,std::back_inserter(limbs),64,false); long long count=0; for(auto limb:limbs) count+=__builtin_popcountll(limb); return count; }
    std::size_t hash_value() const { if (!big_) return std::hash<long long>{}(small_); return boost::multiprecision::hash_value(*big_); }
#define INT_MIXED_OPERATOR(OP, RET) \
template<class T> requires std::is_integral_v<T> friend RET operator OP(const Int& a, T b) { return a OP Int(b); } \
template<class T> requires std::is_integral_v<T> friend RET operator OP(T a, const Int& b) { return Int(a) OP b; }
INT_MIXED_OPERATOR(+, Int)
INT_MIXED_OPERATOR(-, Int)
INT_MIXED_OPERATOR(*, Int)
INT_MIXED_OPERATOR(/, Int)
INT_MIXED_OPERATOR(%, Int)
INT_MIXED_OPERATOR(&, Int)
INT_MIXED_OPERATOR(|, Int)
INT_MIXED_OPERATOR(^, Int)
INT_MIXED_OPERATOR(<<, Int)
INT_MIXED_OPERATOR(>>, Int)
INT_MIXED_OPERATOR(==, bool)
INT_MIXED_OPERATOR(!=, bool)
INT_MIXED_OPERATOR(<, bool)
INT_MIXED_OPERATOR(<=, bool)
INT_MIXED_OPERATOR(>, bool)
INT_MIXED_OPERATOR(>=, bool)
#undef INT_MIXED_OPERATOR
#define INT_FLOAT_OPERATOR(OP) \
template<class T> requires std::is_floating_point_v<T> friend auto operator OP(const Int& a, T b) -> decltype(b OP b) { return static_cast<T>(a) OP b; } \
template<class T> requires std::is_floating_point_v<T> friend auto operator OP(T a, const Int& b) -> decltype(a OP a) { return a OP static_cast<T>(b); }
INT_FLOAT_OPERATOR(+)
INT_FLOAT_OPERATOR(-)
INT_FLOAT_OPERATOR(*)
INT_FLOAT_OPERATOR(/)
INT_FLOAT_OPERATOR(==)
INT_FLOAT_OPERATOR(!=)
INT_FLOAT_OPERATOR(<)
INT_FLOAT_OPERATOR(<=)
INT_FLOAT_OPERATOR(>)
INT_FLOAT_OPERATOR(>=)
#undef INT_FLOAT_OPERATOR

};

namespace std {
    template<> struct hash<Int> { size_t operator()(const Int& value) const { return value.hash_value(); } };
    template<class T> requires is_arithmetic_v<T> struct common_type<Int,T> { using type = conditional_t<is_floating_point_v<T>,T,Int>; };
    template<class T> requires is_arithmetic_v<T> struct common_type<T,Int> { using type = conditional_t<is_floating_point_v<T>,T,Int>; };
}
inline Int abs(const Int& value) { return value < 0 ? -value : value; }
inline std::string to_string(const Int& value) { return value.str(); }
'''
