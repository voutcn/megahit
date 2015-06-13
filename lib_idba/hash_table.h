/**
 * @file hash_table.h
 * @brief HashTableST Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-03
 */

#ifndef __LIB_IDBA_HASH_TABLE_H_

#define __LIB_IDBA_HASH_TABLE_H_

#include <stdint.h>

#include <cstddef>
#include <functional>
#include <istream>
#include <ostream>
#include <stdexcept>

#include "lib_idba/pool.h"
#include "hash.h"
#include "functional.h"


template <typename T>
struct HashTableSTNode
{
    HashTableSTNode<T> *next;
    T value;
};

template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
class HashTableST;

template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
class HashTableSTIterator;

template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
class HashTableSTConstIterator;

template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
class HashTableSTIterator
{
public:
    typedef Key key_type;
    typedef Value value_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef HashTableSTNode<value_type> node_type;
    typedef HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey> hash_table_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef HashTableSTIterator<Value, Key, HashFunc, ExtractKey, EqualKey> iterator;

    HashTableSTIterator(const hash_table_type *owner = NULL, node_type *current = NULL)
        : owner_(owner), current_(current) {}
    HashTableSTIterator(const iterator &iter)
        : owner_(iter.owner_), current_(iter.current_) {}

    const iterator &operator =(const iterator &iter)
    { owner_ = iter.owner_; current_ = iter.current_; return *this; }

    bool operator ==(const iterator &iter) const
    { return current_ == iter.current_; }
    bool operator !=(const iterator &iter) const
    { return current_ != iter.current_; }

    reference operator*() const { return current_->value; }
    pointer operator->() const { return &current_->value; }

    const iterator &operator ++()
    { increment();  return *this; }
    iterator operator ++(int)
    { iterator tmp(*this); increment(); return tmp; }

private:
    void increment()
    {
        if (current_ != NULL)
        {
            if (current_->next)
                current_ = current_->next;
            else
            {
                uint64_t index = owner_->bucket_index_value(current_->value);
                current_ = current_->next;
                while (current_ == NULL && ++index < owner_->bucket_count())
                    current_ = owner_->buckets_[index];
            }
        }
    }

    const hash_table_type *owner_;
    node_type *current_;
};

template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
class HashTableSTConstIterator
{
public:
    typedef Key key_type;
    typedef Value value_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef HashTableSTNode<value_type> node_type;
    typedef HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey> hash_table_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef HashTableSTConstIterator<Value, Key, HashFunc, ExtractKey, EqualKey> const_iterator;

    HashTableSTConstIterator(const hash_table_type *owner = NULL, const node_type *current = NULL)
        : owner_(owner), current_(current) {}
    HashTableSTConstIterator(const const_iterator &iter)
        : owner_(iter.owner_), current_(iter.current_) {}

    const const_iterator &operator =(const const_iterator &iter)
    { owner_ = iter.owner_; current_ = iter.current_; return *this; }

    bool operator ==(const const_iterator &iter) const
    { return current_ == iter.current_; }
    bool operator !=(const const_iterator &iter) const
    { return current_ != iter.current_; }

    const_reference operator*() const { return current_->value; }
    const_pointer operator->() const { return &current_->value; }

    const const_iterator &operator ++()
    { increment();  return *this; }
    const_iterator operator ++(int)
    { const_iterator tmp(*this); increment(); return tmp; }

private:
    void increment()
    {
        if (current_ != NULL)
        {
            if (current_->next)
                current_ = current_->next;
            else
            {
                uint64_t index = owner_->bucket_index_value(current_->value);
                current_ = current_->next;
                while (current_ == NULL && ++index < owner_->bucket_count())
                    current_ = owner_->buckets_[index];
            }
        }
    }

    const hash_table_type *owner_;
    const node_type *current_;
};


/**
 * @brief It is parallel hash table. All insertion/delection operations can be
 * done in parallel. The table size grows automatically, if the number elements
 * exceed the twice of the number of buckets. 
 *
 * @tparam Value
 * @tparam Key
 * @tparam HashFunc
 */
template <typename Value, typename Key, typename HashFunc = Hash<Key>, 
         typename ExtractKey = GetKey<Key, Value>, typename EqualKey = std::equal_to<Key> >
class HashTableST
{
public:
    typedef Key key_type;
    typedef Value value_type;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;

    typedef HashFunc hash_func_type;
    typedef ExtractKey get_key_func_type;
    typedef EqualKey key_equal_func_type;

    typedef HashTableSTNode<value_type> node_type;
    typedef HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey> hash_table_type;
    typedef HashTableSTIterator<Value, Key, HashFunc, ExtractKey, EqualKey> iterator;
    typedef HashTableSTConstIterator<Value, Key, HashFunc, ExtractKey, EqualKey> const_iterator;
    typedef PoolST<node_type> pool_type;

    friend class HashTableSTIterator<Value, Key, HashFunc, ExtractKey, EqualKey>;
    friend class HashTableSTConstIterator<Value, Key, HashFunc, ExtractKey, EqualKey>;

    template <typename Value_, typename Key_, typename HashFunc_,
             typename ExtractKey_, typename EqualKey_>
    friend std::ostream &operator <<(std::ostream &os, 
            HashTableST<Value_, Key_, HashFunc_, ExtractKey_, EqualKey_> &hash_table);

    template <typename Value_, typename Key_, typename HashFunc_,
             typename ExtractKey_, typename EqualKey_>
    friend std::istream &operator >>(std::istream &os, 
            HashTableST<Value_, Key_, HashFunc_, ExtractKey_, EqualKey_> &hash_table);

    static const uint64_t kNumBucketLocks = (1 << 12);
    static const uint64_t kDefaultNumBuckets = (1 << 12);

    explicit HashTableST(const hash_func_type &hash = hash_func_type(),
            const get_key_func_type &get_key = get_key_func_type(),
            const key_equal_func_type &key_equal = key_equal_func_type())
        : hash_(hash), get_key_(get_key), key_equal_(key_equal)
    { 
        size_ = 0;
        rehash(kDefaultNumBuckets); 
    }

    HashTableST(const hash_table_type &hash_table)
        : hash_(hash_table.hash_),
        get_key_(hash_table.get_key_),
        key_equal_(hash_table.key_equal_)
    {
        size_ = 0;
        assign(hash_table);
    }

    ~HashTableST()
    {
        clear();
    }

    const hash_table_type &operator =(const hash_table_type &hash_table)
    { return assign(hash_table); }

    const hash_table_type &assign(const hash_table_type &hash_table)
    {
        if (this == &hash_table)
            return *this;

        clear();
        rehash(hash_table.buckets_.size());

        for (int64_t i = 0; i < (int64_t)hash_table.buckets_.size(); ++i)
        {
            node_type *prev = NULL;
            for (node_type *node = hash_table.buckets_[i]; node; node = node->next)
            {
                node_type *p = pool_.construct();
                p->value = node->value;
                p->next = NULL;

                if (prev == NULL)
                    buckets_[i] = p;
                else
                    prev->next = p;

                prev = p;
            }
        }

        return *this;
    }

    iterator begin()
    {
        for (unsigned i = 0; i < buckets_.size(); ++i)
        {
            if (buckets_[i])
                return iterator(this, buckets_[i]);
        }
        return iterator();
    } 

    const_iterator begin() const
    {
        for (unsigned i = 0; i < buckets_.size(); ++i)
        {
            if (buckets_[i])
                return const_iterator(this, buckets_[i]);
        }
        return const_iterator();
    }

    iterator end()
    { return iterator(); }

    const_iterator end() const
    { return const_iterator(); }

    std::pair<iterator, bool> insert_unique(const value_type &value)
    {
        rehash_if_needed(size_);

        uint64_t hash_value = hash(value);
        lock_bucket(hash_value);
        uint64_t index = bucket_index(hash_value);

        for (node_type *node = buckets_[index]; node; node = node->next)
        {
            if (key_equal_(get_key_(node->value), get_key_(value)))
            {
                unlock_bucket(hash_value);
                return std::pair<iterator, bool>(iterator(this, node), false);
            }
        }

        node_type *p= pool_.construct();
        p->value = value;
        p->next = buckets_[index];
        buckets_[index] = p;
        ++size_;
        unlock_bucket(hash_value);

        return std::pair<iterator, bool>(iterator(this, p), true);
    }

    iterator find(const key_type &key)
    {
        uint64_t hash_value = hash_key(key);
        lock_bucket(hash_value);
        uint64_t index = bucket_index_key(key);
        for (node_type *node = buckets_[index]; node; node = node->next)
        {
            if (key_equal_(key, get_key_(node->value)))
            {
                unlock_bucket(hash_value);
                return iterator(this, node);
            }
        }
        unlock_bucket(hash_value);
        return iterator();
    }

    const_iterator find(const key_type &key) const
    {
        uint64_t index = bucket_index_key(key);
        for (node_type *node = buckets_[index]; node; node = node->next)
        {
            if (key_equal_(key, get_key_(node->value)))
            {
                return const_iterator(this, node);
            }
        }
        return const_iterator();
    }

    reference find_or_insert(const value_type &value)
    {
        rehash_if_needed(size_);

        uint64_t hash_value = hash(value);
        lock_bucket(hash_value);
        uint64_t index = bucket_index(hash_value);

        for (node_type *node = buckets_[index]; node; node = node->next)
        {
            if (key_equal_(get_key_(node->value), get_key_(value)))
            {
                unlock_bucket(hash_value);
                return node->value;
            }
        }

        node_type *p= pool_.construct();
        p->value = value;
        p->next = buckets_[index];
        buckets_[index] = p;
        ++size_;
        unlock_bucket(hash_value);

        return p->value;
    }

    size_type remove(const key_type &key)
    {
        uint64_t num_removed_nodes = 0;

        uint64_t hash_value = hash_key(key);
        lock_bucket(hash_value);
        uint64_t index = bucket_index(hash_value);

        node_type *prev = NULL;
        node_type *node = buckets_[index];
        while (node)
        {
            if (key_equal_(key, get_key_(node->value)))
            {
                if (prev == NULL)
                    buckets_[index] = node->next;
                else
                    prev->next = node->next;

                node_type *p = node;
                node = node->next;
                pool_.destroy(p);

                ++num_removed_nodes;
            }
            else
            {
                prev = node;
                node = node->next;
            }
        }
        unlock_bucket(hash_value);

        size_ -= num_removed_nodes;

        return num_removed_nodes;
    }

    template <typename Predicator>
    size_type remove_if(const Predicator &predicator)
    {
        uint64_t num_removed_nodes = 0;

        for (int64_t index = 0; index < (int64_t)buckets_.size(); ++index)
        {
            lock_bucket(index);

            node_type *prev = NULL;
            node_type *node = buckets_[index];
            while (node)
            {
                if (predicator(node->value))
                {
                    if (prev == NULL)
                        buckets_[index] = node->next;
                    else
                        prev->next = node->next;

                    node_type *p = node;
                    node = node->next;
                    pool_.destroy(p);

                    ++num_removed_nodes;
                }
                else
                {
                    prev = node;
                    node = node->next;
                }
            }
            unlock_bucket(index);
        }

        size_ -= num_removed_nodes;

        return num_removed_nodes;
    }

    template <typename UnaryProc>
    UnaryProc &for_each(UnaryProc &op)
    {
        for (int64_t i = 0; i < (int64_t)buckets_.size(); ++i)
        {
            for (node_type *node = buckets_[i]; node; node = node->next)
                op(node->value);
        }
        return op;
    }

    template <typename UnaryProc>
    UnaryProc &for_each(UnaryProc &op) const
    {
        for (int64_t i = 0; i < (int64_t)buckets_.size(); ++i)
        {
            for (node_type *node = buckets_[i]; node; node = node->next)
                op(node->value);
        }
        return op;
    }

    uint64_t hash(const value_type &value) const
    { return hash_(get_key_(value)); }

    uint64_t bucket_index(uint64_t hash_value) const
    { return hash_value & (buckets_.size() -1); }

    uint64_t hash_key(const key_type &key) const
    { return hash_(key); }

    uint64_t bucket_index_value(const value_type &value) const
    { return hash_(get_key_(value)) & (buckets_.size() - 1); }

    uint64_t bucket_index_key(const key_type &key) const
    { return hash_(key) & (buckets_.size() - 1); }

    const hash_func_type &hash_func() const
    { return hash_; }
    const get_key_func_type &get_key_func() const
    { return get_key_; }
    const key_equal_func_type &key_equal_func() const
    { return key_equal_; }

    size_type bucket_count() const { return buckets_.size(); }

    void reserve(size_type capacity)
    { rehash_if_needed(capacity); }

    void swap(hash_table_type &hash_table)
    {
        if (this != &hash_table)
        {
            std::swap(hash_, hash_table.hash_);
            std::swap(get_key_, hash_table.get_key_);
            std::swap(key_equal_, hash_table.key_equal_);

            pool_.swap(hash_table.pool_);
            buckets_.swap(hash_table.buckets_);
            std::swap(size_, hash_table.size_);
        }
    }

    size_type size() const { return size_; }
    bool empty() const { return size_ == 0; }

    void clear()
    {
        size_ = 0;

        for (int64_t i = 0; i < (int64_t)buckets_.size(); ++i)
        {
            node_type *node = buckets_[i];
            while (node)
            {
                node_type *p = node;
                node = node->next;
                pool_.destroy(p);
            }
            buckets_[i] = NULL;
        }
        pool_.clear();
    }

private:
    void lock_bucket(uint64_t hash_value)
    { }
    void unlock_bucket(uint64_t hash_value)
    { }

    void rehash_if_needed(size_type capacity)
    {
        if (capacity > buckets_.size() * 2)
        {
            size_type new_num_buckets = buckets_.size();
            while (capacity > new_num_buckets * 2)
                new_num_buckets *= 2;
            rehash(new_num_buckets);
        }
    }

    void rehash(uint64_t new_num_buckets)
    {
        if ((new_num_buckets & (new_num_buckets-1)) != 0)
            throw std::logic_error("HashTableST::rehash() invalid number of buckets");

        if (new_num_buckets == buckets_.size())
            return;

        std::vector<node_type *> old_buckets(new_num_buckets, NULL);
        old_buckets.swap(buckets_);

        if (new_num_buckets > old_buckets.size())
        {
            for (int64_t i = 0; i < (int64_t)old_buckets.size(); ++i)
            {
                node_type *node = old_buckets[i];
                while (node)
                {
                    node_type *next = node->next;
                    uint64_t index = bucket_index_value(node->value);
                    node->next = buckets_[index];
                    buckets_[index] = node;
                    node = next;
                }
            }
        }
        else
        {
            for (int64_t i = 0; i < (int64_t)old_buckets.size(); ++i)
            {
                node_type *node = old_buckets[i];
                while (node)
                {
                    uint64_t hash_value = hash(node->value);
                    //lock_bucket(hash_value);
                    uint64_t index = bucket_index(hash_value);
                    node_type *next = node->next;
                    node->next = buckets_[index];
                    buckets_[index] = node;
                    node = next;
                    //unlock_bucket(hash_value);
                }
            }
        }
    }

    hash_func_type hash_;
    get_key_func_type get_key_;
    key_equal_func_type key_equal_;

    PoolST<node_type> pool_;
    std::vector<node_type *> buckets_;
    uint64_t size_;
};

template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
std::istream &operator >>(std::istream &is, 
        HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey> &hash_table)
{
    hash_table.clear();
    Value value;
    while (is.read((char *)&value, sizeof(Value)))
        hash_table.insert_unique(value);
    return is;
}

template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
std::ostream &operator <<(std::ostream &os, 
        HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey> &hash_table)
{
    typename HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey>::iterator iter;
    for (iter = hash_table.begin(); iter != hash_table.end(); ++iter)
    {
        os.write((char *)&*iter, sizeof(Value));
    }
    return os;
}

namespace std
{
template <typename Value, typename Key, typename HashFunc,
         typename ExtractKey, typename EqualKey>
inline void swap(HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey> &x,
     HashTableST<Value, Key, HashFunc, ExtractKey, EqualKey> &y)
{ x.swap(y); }
}

#endif
