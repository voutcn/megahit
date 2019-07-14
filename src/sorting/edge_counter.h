//
// Created by dinghua.li on 6/30/19.
//

#ifndef MEGAHIT_EDGE_COUNTER_H
#define MEGAHIT_EDGE_COUNTER_H

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include "sdbg/sdbg_def.h"

class EdgeMultiplicityRecorder {
 public:
  EdgeMultiplicityRecorder() = default;

  void SetNumThreads(unsigned n) {
    counters_.resize(n);
    for (auto &counter : counters_) {
      std::fill(counter.begin(), counter.end(), 0);
    }
  }

  size_t size_in_byte() const {
    return (kMaxMul + 1) * counters_.size() * sizeof(int64_t);
  }

  template <typename T>
  void Add(T multiplicity, unsigned thread_id) {
    ++counters_[thread_id][std::min(static_cast<T>(kMaxMul), multiplicity)];
  }

  int64_t GetNumSolidEdges(int solid_threshold) const {
    int64_t sum = 0;
    for (const auto &counter : counters_) {
      for (int i = solid_threshold; i <= kMaxMul; ++i) {
        sum += counter[i];
      }
    }
    return sum;
  }

  void DumpStat(std::ostream &os) const {
    for (int i = 1; i <= kMaxMul; ++i) {
      int64_t sum = 0;
      for (const auto &counter : counters_) {
        sum += counter[i];
      }
      os << i << ' ' << sum << '\n';
    }
  }

 private:
  std::vector<std::array<int64_t, kMaxMul + 1>> counters_;
};

#endif  // MEGAHIT_EDGE_COUNTER_H
