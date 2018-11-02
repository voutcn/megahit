/**
 * @file pool.h
 * @brief Pool Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-23
 */

#ifndef __LIB_IDBA_POOL_H_

#define __LIB_IDBA_POOL_H_

#include <stdint.h>

#include <deque>
#include <vector>

template <typename T>
struct ChunkST {
    ChunkST(T *address = NULL, uint32_t size = 0) {
        this->address = address;
        this->size = size;
    }

    T *address;
    uint32_t size;
};

template <typename T>
struct BufferST {
    BufferST(T *address = NULL, uint32_t size = 0, uint32_t index = 0) {
        this->address = address;
        this->size = size;
        this->index = index;
    }

    T *address;
    uint32_t size;
    uint32_t index;
};

template <typename T, typename Allocator = std::allocator<T> >
class PoolST {
  public:
    typedef T value_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef Allocator allocator_type;
    typedef ChunkST<T> chunk_type;
    typedef BufferST<T> buffer_type;
    typedef PoolST<T> pool_type;

//    static const uint32_t kMaxChunkSTSize = (1 << 12);
//    static const uint32_t kMinChunkSTSize = (1 << 12);
    static const uint32_t kMaxChunkSTSize = (1 << 20);
    static const uint32_t kMinChunkSTSize = (1 << 8);


    PoolST() {
        heads_.resize(1, (pointer)0);
        buffers_.resize(1, buffer_type());
        chunk_size_ = kMinChunkSTSize;
    }
    ~PoolST() {
        clear();
    }

    pointer allocate() {
        int thread_id = 0;
        if (heads_[thread_id] != NULL) {
            pointer p = heads_[thread_id];
            heads_[thread_id] = *(pointer *)heads_[thread_id];
            return p;
        } else {
            buffer_type &buffer = buffers_[thread_id];
            if (buffer.index == buffer.size) {
                uint32_t size = chunk_size_;
                if (chunk_size_ < kMaxChunkSTSize)
                    chunk_size_ <<= 1;

                pointer p = alloc_.allocate(size);

                chunks_.push_back(chunk_type(p, size));

                buffer.address = p;
                buffer.size = size;
                buffer.index = 0;
            }

            return buffer.address + buffer.index++;
        }
    }

    void deallocate(pointer p) {
        int thread_id = 0;
        *(pointer *)p = heads_[thread_id];
        heads_[thread_id] = p;
    }

    pointer construct() {
        pointer p = allocate();
        new ((void *)p)value_type();
        return p;
    }

    pointer construct(const_reference x) {
        pointer p = allocate();
        new ((void *)p)value_type(x);
        return p;
    }

    void destroy(pointer p) {
        ((value_type*)p)->~value_type();
    }

    void swap(PoolST<value_type> &pool) {
        if (this != &pool) {
            heads_.swap(pool.heads_);
            buffers_.swap(pool.buffers_);
            chunks_.swap(pool.chunks_);
            std::swap(chunk_size_, pool.chunk_size_);
            std::swap(alloc_, pool.alloc_);
        }
    }

    void clear() {
        for (unsigned i = 0; i < chunks_.size(); ++i)
            alloc_.deallocate(chunks_[i].address, chunks_[i].size);
        chunks_.resize(0);
        fill(heads_.begin(), heads_.end(), (pointer)0);
        fill(buffers_.begin(), buffers_.end(), buffer_type());
        chunk_size_ = kMinChunkSTSize;
    }

  private:
    PoolST(const pool_type &);
    const pool_type &operator =(const pool_type &);

    std::vector<pointer> heads_;
    std::vector<buffer_type> buffers_;
    std::deque<chunk_type> chunks_;
    uint32_t chunk_size_;
    allocator_type alloc_;
};

namespace std {
template <typename T> inline void swap(PoolST<T> &x, PoolST<T> &y) {
    x.swap(y);
}
}

#endif

