//
// Created by vout on 5/25/19.
//

#ifndef MEGAHIT_SPINLOCK_H
#define MEGAHIT_SPINLOCK_H

#include <atomic>

class SpinLock {
 public:
  void lock() {
    while (flag_.test_and_set(std::memory_order_acquire)) {}
  }
  void unlock() {
    flag_.clear(std::memory_order_release);
  }

 private:
  std::atomic_flag flag_ = ATOMIC_FLAG_INIT;
};
#endif //MEGAHIT_SPINLOCK_H
