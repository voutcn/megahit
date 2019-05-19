//
// Created by vout on 3/7/19.
//

#ifndef MEGAHIT_ATOMIC_WRAPPER_H
#define MEGAHIT_ATOMIC_WRAPPER_H

#include <atomic>

/*!
 * @brief a wrapper for std::Atomic. std::Atomic do not support copy and move
 * constructor, so this wrapper is used to make suitable to std::vector
 * @tparam T the underlying type of the atomic struct
 */
template <typename T>
struct AtomicWrapper {
  std::atomic<T> v;
  AtomicWrapper(T a = T()) : v(a) {}
  AtomicWrapper(const AtomicWrapper &rhs) : v(rhs.v.load()) {}
  AtomicWrapper &operator=(const AtomicWrapper &rhs) {
    v.store(rhs.v.load());
    return *this;
  }
};

#endif  // MEGAHIT_ATOMIC_WRAPPER_H
