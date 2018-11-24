//
// Created by vout on 11/23/18.
//

#ifndef MEGAHIT_JUNCION_INDEX_H
#define MEGAHIT_JUNCION_INDEX_H

#include <mutex>
#include <sdbg/sdbg_def.h>
#include <utils.h>
#include "sequence/kmer_plus.h"
#include "sparsepp/spp.h"
#include "sequence_package.h"

template<class KmerType>
class JunctionIndex {
 public:
  struct JunctionInfo {
    uint64_t ext_seq : 58;
    unsigned ext_len : 6;
    float mul;
  } __attribute__((packed));
  using Junction = KmerPlus<KmerType::kNumWords, typename KmerType::word_type, JunctionInfo>;
  using hash_set = spp::sparse_hash_set<Junction, KmerHash>;
 public:
  JunctionIndex(unsigned k, unsigned step) : k_(k), step_(step) {}
  size_t size() const { return hash_index_.size(); }

  void FeedBatchContigs(SequencePackage &seq_pkg, const std::vector<float> &mul) {
    std::mutex lock;
#pragma omp parallel for
    for (size_t i = 0; i < seq_pkg.size(); ++i) {
      size_t seq_len = seq_pkg.length(i);
      if (seq_len < k_ + 1) {
        continue;
      }
      for (int strand = 0; strand < 2; ++strand) {
        auto get_jth_char = [&seq_pkg, i, strand, seq_len](unsigned j) -> uint8_t {
          uint8_t c = seq_pkg.get_base(i, strand == 0 ? j : (seq_len - 1 - j));
          return strand == 0 ? c : 3 ^ c;
        };

        Junction jct;
        for (unsigned j = 0; j < k_ + 1; ++j) {
          jct.kmer.ShiftAppend(get_jth_char(j), k_ + 1);
        }
        if (jct.kmer.IsPalindrome(k_ + 1)) {
          continue;
        }

        unsigned ext_len = std::min(static_cast<size_t>(step_ - 1), seq_len - (k_ + 1));
        uint64_t ext_seq = 0;
        for (unsigned j = 0; j < ext_len; ++j) {
          ext_seq |= uint64_t(get_jth_char(k_ + 1 + j)) << j * 2;
        }
        jct.aux.ext_len = ext_len;
        jct.aux.ext_seq = ext_seq;
        jct.aux.mul = mul[i];

        {
          std::lock_guard<std::mutex> lk(lock);
          auto iter = hash_index_.find(jct);
          if (iter != hash_index_.end()) {
            if (iter->aux.ext_len < ext_len || (iter->aux.ext_len == ext_len && iter->aux.ext_seq < ext_seq)) {
              hash_index_.erase(iter);
              hash_index_.insert(jct);
            }
          } else {
            hash_index_.insert(jct);
          }
        }
        if (seq_len == k_ + 1) {
          break;
        }
      }
    }
  }

  template<class CollectorType>
  size_t FindNextKmersFromRead(
      SequencePackage &seq_pkg, size_t seq_id,
      CollectorType *out) {
    size_t length = seq_pkg.length(seq_id);
    if (length < k_ + step_ + 1) {
      return 0;
    }

    size_t num_success = 0;
    std::vector<bool> kmer_exist(length, false);
    std::vector<float> kmer_mul(length, 0);

    Junction jct, rjct;
    auto &kmer = jct.kmer;
    auto &rkmer = rjct.kmer;

    for (unsigned j = 0; j < k_ + 1; ++j) {
      kmer.ShiftAppend(seq_pkg.get_base(seq_id, j), k_ + 1);
    }
    rkmer = kmer;
    rkmer.ReverseComplement(k_ + 1);

    unsigned cur_pos = 0;
    while (cur_pos + k_ + 1 <= length) {
      unsigned next_pos = cur_pos + 1;

      if (!kmer_exist[cur_pos]) {
        auto iter = hash_index_.find(jct);
        if (iter != hash_index_.end()) {
          kmer_exist[cur_pos] = true;
          uint64_t ext_seq = iter->aux.ext_seq;
          unsigned ext_len = iter->aux.ext_len;
          float mul = iter->aux.mul;
          kmer_mul[cur_pos] = mul;

          for (unsigned j = 0; j < ext_len && cur_pos + k_ + 1 + j < length; ++j, ++next_pos) {
            if (seq_pkg.get_base(seq_id, cur_pos + k_ + 1 + j) == ((ext_seq >> j * 2) & 3)) {
              kmer_exist[cur_pos + j + 1] = true;
              kmer_mul[cur_pos + j + 1] = mul;
            } else {
              break;
            }
          }
        } else if ((iter = hash_index_.find(rjct)) != hash_index_.end()) {
          kmer_exist[cur_pos] = true;
          uint64_t ext_seq = iter->aux.ext_seq;
          unsigned ext_len = iter->aux.ext_len;
          float mul = iter->aux.mul;
          kmer_mul[cur_pos] = mul;

          for (unsigned j = 0; j < ext_len && cur_pos >= j + 1; ++j) {
            if ((3 ^ seq_pkg.get_base(seq_id, cur_pos - 1 - j)) == ((ext_seq >> j * 2) & 3)) {
              kmer_exist[cur_pos - 1 - j] = true;
              kmer_mul[cur_pos - 1 - j] = mul;
            } else {
              break;
            }
          }
        }
      }

      if (next_pos + k_ + 1 <= length) {
        while (cur_pos < next_pos) {
          ++cur_pos;
          uint8_t c = seq_pkg.get_base(seq_id, cur_pos + k_);
          kmer.ShiftAppend(c, k_ + 1);
          rkmer.ShiftPreappend(3 ^ c, k_ + 1);
        }
      } else {
        break;
      }
    }

    for (int j = 1; j + k_ + 1 <= length; ++j) {
      kmer_mul[j] += kmer_mul[j - 1];
    }

    typename CollectorType::kmer_type new_kmer, new_rkmer;

    for (unsigned accumulated_len = 0, j = 0, end_pos = 0; j + k_ < length; ++j) {
      accumulated_len = kmer_exist[j] ? accumulated_len + 1 : 0;
      if (accumulated_len >= step_ + 1) {
        unsigned target_end = j + k_ + 1;
        if (end_pos + 8 < target_end) {
          while (end_pos < target_end) {
            auto c = seq_pkg.get_base(seq_id, end_pos);
            new_kmer.ShiftAppend(c, k_ + step_ + 1);
            new_rkmer.ShiftPreappend(3 ^ c, k_ + step_ + 1);
            end_pos++;
          }
        } else {
          if (end_pos + k_ + step_ + 1 < target_end) {
            end_pos = target_end - (k_ + step_ + 1);
          }
          while (end_pos < target_end) {
            auto c = seq_pkg.get_base(seq_id, end_pos);
            new_kmer.ShiftAppend(c, k_ + step_ + 1);
            end_pos++;
          }
          new_rkmer = new_kmer;
          new_rkmer.ReverseComplement(k_ + step_ + 1);
        }
        float mul = (kmer_mul[j] - (j >= step_ + 1 ? kmer_mul[j - (step_ + 1)] : 0)) / (step_ + 1);
        assert(mul <= kMaxMul + 1);
        out->insert(new_kmer < new_rkmer ? new_kmer : new_rkmer,
                    static_cast<mul_t>(std::min(kMaxMul, static_cast<int>(mul + 0.5))));
        num_success++;
      }
    }
    return num_success;
  }
 private:
  hash_set hash_index_;
  unsigned k_{};
  unsigned step_{};
};

#endif //MEGAHIT_JUNCION_INDEX_H
