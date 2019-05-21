//
// Created by vout on 5/20/19.
//

#ifndef MEGAHIT_HASHTABLE_WRAPPER_H
#define MEGAHIT_HASHTABLE_WRAPPER_H

#include "idba/hash.h"
#include "sparsepp/spp.h"
#include "idba/functional.h"
#include "idba/khash.h"

template<typename Value, typename Key>
class SppHashTable :
    public spp::sparse_hashtable<Value,
                                 Key,
                                 Hash<Key>,
                                 GetKey<Key, Value>,
                                 SetKey<Value>,
                                 std::equal_to<Key>,
                                 SPP_DEFAULT_ALLOCATOR<Value>> {
 public:
  Value &find_or_insert(const Value &val) {
    auto ret = this->insert(val);
    return *ret.first;
  }

  template<typename UnaryProc>
  void for_each(UnaryProc &op) {
    for (auto i = this->begin(); i != this->end();
         ++i) {
      op(*i);
    }
  }

  template<typename UnaryProc>
  void for_each(UnaryProc &op) const {
    for (auto i = this->begin(); i != this->end();
         ++i) {
      op(*i);
    }
  }

  void reserve(size_t sz) {}
};


template<typename Value, typename Key>
class KHashTable {
 public:
  struct HashViaValue {
    uint64_t operator() (const Value &value) {
      return GetKey<Key, Value>()(value).hash();
    }
  };

  struct KeyEqual {
    uint64_t operator() (const Value &val1, const Value &val2) {
      return GetKey<Key, Value>()(val1) == GetKey<Key, Value>()(val2);
    }
  };

  KHashTable() : h_(new klib::KHash<Value, HashViaValue, KeyEqual, uint64_t>()) {}

  Value &find_or_insert(const Value &val) {
    int ret;
    return h_->at(h_->put(val, &ret));
  }

  template<typename UnaryProc>
  void for_each(UnaryProc &op) {
    for (auto i = h_->begin(); i != h_->end(); ++i) {
      if (h_->exist(i)) {
        op(h_->at(i));
      }
    }
  }

  void reserve(size_t sz) { h_->resize(sz); }

  Value* find(const Value &val) {
    auto at = h_->get(val);
    if (at == h_->end()) {
      return nullptr;
    }
    return &(h_->at(at));
  }

  Value* end() const {
    return nullptr;
  }

  void swap(KHashTable &kh) {
    std::swap(h_, kh.h_);
  }
  void clear() {
    h_.reset(new klib::KHash<Value, HashViaValue, KeyEqual, uint64_t>());
  }
  size_t size() const {
    return h_->size();
  }

 private:
  std::unique_ptr<klib::KHash<Value, HashViaValue, KeyEqual, uint64_t>> h_;
};

#endif //MEGAHIT_HASHTABLE_WRAPPER_H
