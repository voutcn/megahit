//
// Created by vout on 5/20/19.
//

#ifndef MEGAHIT_HASHTABLE_WRAPPER_H
#define MEGAHIT_HASHTABLE_WRAPPER_H

#include "idba/hash.h"
#include "sparsepp/spp.h"
#include "idba/functional.h"

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
    for (auto i = this->begin(); i != this->end(); ++i) {
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

#endif //MEGAHIT_HASHTABLE_WRAPPER_H
