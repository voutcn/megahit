/**
 * @file atomic_integer.h
 * @brief AtomicInteger Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-05
 */

#ifndef __BASIC_ATOMIC_INTEGER_H_

#define __BASIC_ATOMIC_INTEGER_H_

#include <algorithm>

/**
 * @brief It is a integer-like class which is used for atomic access.
 *
 * @tparam T
 */
template <typename T>
class AtomicInteger
{
public:
    AtomicInteger(): value_(0) {}
    explicit AtomicInteger(T x): value_(x) {}
    AtomicInteger(const AtomicInteger<T> &x): value_(x.value_) {}

    operator T() const { return value_; }

    T operator =(const AtomicInteger<T> &x) { value_ = x.value_; return *this; }
    T operator =(T x) { value_ = x; return *this; }

    bool operator <(const AtomicInteger<T> &x) const { return value_ < x.value_; }
    bool operator >(const AtomicInteger<T> &x) const { return value_ > x.value_; }
    bool operator ==(const AtomicInteger<T> &x) const { return value_ == x.value_; }
    bool operator !=(const AtomicInteger<T> &x) const { return value_ != x.value_; }

    T operator += (T x) { return __sync_add_and_fetch(&value_, x); }
    T operator -= (T x) { return __sync_sub_and_fetch(&value_, x); }
    T operator |= (T x) { return __sync_or_and_fetch(&value_, x); }
    T operator &= (T x) { return __sync_and_and_fetch(&value_, x); }
    T operator ^= (T x) { return __sync_xor_and_fetch(&value_, x); }

    T operator ++() { return __sync_add_and_fetch(&value_, 1); }
    T operator ++(int) { return __sync_fetch_and_add(&value_, 1); }
    T operator --() { return __sync_sub_and_fetch(&value_, 1); }
    T operator --(int) { return __sync_fetch_and_sub(&value_, 1); }

    bool CompareAndSet(T old_value, T new_value)
    { return __sync_bool_compare_and_swap(&value_, old_value, new_value); }

    void swap(AtomicInteger &x) 
    { if (this != &x) std::swap(value_, x.value_); }

private:
    T value_;
};

namespace std
{
template <typename T> inline void swap(AtomicInteger<T> &x, AtomicInteger<T> &y) { x.swap(y); }
}

#endif

