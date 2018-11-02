/**
 * @file hash_map.h
 * @brief HashMap Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-24
 * @last modified by Dinghua Li
 * @date 2014-10-06
 */

#ifndef __CONTAINER_HASH_MAP_H_

#define __CONTAINER_HASH_MAP_H_

#include "functional.h"
#include "hash.h"
#include "hash_table.h"

#include <functional>


/**
 * @brief It is a parallel hash map which has similar inserface as stl map.
 * It is implemented based on parallel hash table (HashTable).
 *
 * @tparam Key
 * @tparam Value
 * @tparam HashFunc
 */
template <typename Key, typename Value, typename HashFunc = Hash<Key>,
          typename EqualKey = std::equal_to<Key> >
class HashMap {
  public:
    typedef HashTable<std::pair<Key, Value>, Key, HashFunc,
            Select1st<std::pair<Key, Value> >, EqualKey> hash_table_type;
    typedef HashMap<Key, Value, HashFunc, EqualKey> hash_map_type;

    typedef typename hash_table_type::key_type key_type;
    typedef typename hash_table_type::value_type value_type;
    typedef typename hash_table_type::size_type size_type;
    typedef typename hash_table_type::difference_type difference_type;

    typedef typename hash_table_type::reference reference;
    typedef typename hash_table_type::const_reference const_reference;
    typedef typename hash_table_type::pointer pointer;
    typedef typename hash_table_type::const_pointer const_pointer;

    typedef Value data_type;

    typedef typename hash_table_type::hash_func_type hash_func_type;
    typedef typename hash_table_type::get_key_func_type get_key_func_type;
    typedef typename hash_table_type::key_equal_func_type key_equal_func_type;

    typedef typename hash_table_type::iterator iterator;
    typedef typename hash_table_type::const_iterator const_iterator;

    explicit HashMap(const hash_func_type &hash = hash_func_type(),
                     const key_equal_func_type key_equal = key_equal_func_type())
        : hash_table_(hash, Select1st<std::pair<Key, Value> >(), key_equal) {
    }

    HashMap(const hash_map_type &hash_map)
        : hash_table_(hash_map.hash_table_) {
    }

    const hash_map_type &operator = (const hash_map_type &hash_map) {
        hash_table_ = hash_map.hash_table_;
        return *this;
    }

    iterator begin() {
        return hash_table_.begin();
    }
    const_iterator begin() const {
        return hash_table_.begin();
    }
    iterator end() {
        return hash_table_.end();
    }
    const_iterator end() const {
        return hash_table_.end();
    }

    std::pair<iterator, bool> insert(const value_type &value) {
        return hash_table_.insert_unique(value);
    }

    iterator find(const key_type &key) {
        return hash_table_.find(key);
    }

    const_iterator find(const value_type &value) const {
        return hash_table_.find(value);
    }

    data_type &operator [](const key_type &key) {
        return hash_table_.find_or_insert(value_type(key, data_type())).second;
    }

    data_type &get_ref_with_lock(const key_type &key) {
        return hash_table_.find_or_insert_with_lock(value_type(key, data_type())).second;
    }

    void unlock(const key_type &key) {
        hash_table_.unlock(value_type(key, data_type()));
    }

    size_type remove(const key_type &key) {
        return hash_table_.remove(key);
    }

    template <typename Predicator>
    size_type remove_if(Predicator &predicator) {
        return hash_table_.remove_if(predicator);
    }

    template <typename UnaryProc>
    UnaryProc &for_each(UnaryProc &op) {
        return hash_table_.for_each(op);
    }

    const hash_func_type &hash_func() const {
        return hash_table_.hash_func();
    }
    const key_equal_func_type &key_equal_func() const {
        return hash_table_.key_equal_func();
    }

    void reserve(size_type capacity) {
        hash_table_.reserve(capacity);
    }

    void swap(hash_map_type &hash_map) {
        if (this != &hash_map) hash_table_.swap(hash_map.hash_table_);
    }

    size_type size() const {
        return hash_table_.size();
    }
    bool empty() const {
        return hash_table_.empty();
    }

    void clear() {
        hash_table_.clear();
    }

  private:
    hash_table_type hash_table_;
};

namespace std {
template <typename Key, typename Value, typename HashFunc, typename EqualKey>
inline void swap(HashMap<Key, Value, HashFunc, EqualKey> &x,
                 HashMap<Key, Value, HashFunc, EqualKey> &y) {
    x.swap(y);
}
}

#endif

