/**
 * @file pool.h
 * @brief Pool Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-23
 */

#ifndef __BASIC_POOL_H_

#define __BASIC_POOL_H_

#include <omp.h>
#include <stdint.h>

#include <deque>
#include <vector>

template <typename T>
struct Chunk {
    Chunk(T *address = NULL, uint32_t size = 0) {
        this->address = address;
        this->size = size;
    }

    T *address;
    uint32_t size;
};

template <typename T>
struct Buffer {
    Buffer(T *address = NULL, uint32_t size = 0, uint32_t index = 0) {
        this->address = address;
        this->size = size;
        this->index = index;
    }

    T *address;
    uint32_t size;
    uint32_t index;
};

template <typename T, typename Allocator = std::allocator<T> >
class Pool {
  public:
    typedef T value_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef Allocator allocator_type;
    typedef Chunk<T> chunk_type;
    typedef Buffer<T> buffer_type;
    typedef Pool<T> pool_type;

    //    static const uint32_t kMaxChunkSize = (1 << 12);
    //    static const uint32_t kMinChunkSize = (1 << 12);
    static const uint32_t kMaxChunkSize = (1 << 20);
    static const uint32_t kMinChunkSize = (1 << 8);


    Pool() {
        omp_init_lock(&lock_alloc_);
        heads_.resize(omp_get_max_threads(), (pointer)0);
        buffers_.resize(omp_get_max_threads(), buffer_type());
        chunk_size_ = kMinChunkSize;
    }
    ~Pool() {
        clear();
        omp_destroy_lock(&lock_alloc_);
    }

    pointer allocate() {
        int thread_id = omp_get_thread_num();

        if (heads_[thread_id] != NULL) {
            pointer p = heads_[thread_id];
            heads_[thread_id] = *(pointer *)heads_[thread_id];
            return p;
        }
        else {
            buffer_type &buffer = buffers_[thread_id];

            if (buffer.index == buffer.size) {
                omp_set_lock(&lock_alloc_);
                uint32_t size = chunk_size_;

                if (chunk_size_ < kMaxChunkSize)
                    chunk_size_ <<= 1;

                omp_unset_lock(&lock_alloc_);

                pointer p = alloc_.allocate(size);

                omp_set_lock(&lock_alloc_);
                chunks_.push_back(chunk_type(p, size));
                // fprintf(stderr, "%p - %p, size: %lu\n", p+1, p, sizeof(*p));
                omp_unset_lock(&lock_alloc_);

                buffer.address = p;
                buffer.size = size;
                buffer.index = 0;
            }

            return buffer.address + buffer.index++;
        }
    }

    void deallocate(pointer p) {
        int thread_id = omp_get_thread_num();
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
        ((value_type *)p)->~value_type();
    }

    void swap(Pool<value_type> &pool) {
        if (this != &pool) {
            heads_.swap(pool.heads_);
            buffers_.swap(pool.buffers_);
            chunks_.swap(pool.chunks_);
            std::swap(chunk_size_, pool.chunk_size_);
            std::swap(alloc_, pool.alloc_);
        }
    }

    void clear() {
        omp_set_lock(&lock_alloc_);

        for (unsigned i = 0; i < chunks_.size(); ++i)
            alloc_.deallocate(chunks_[i].address, chunks_[i].size);

        chunks_.resize(0);
        fill(heads_.begin(), heads_.end(), (pointer)0);
        fill(buffers_.begin(), buffers_.end(), buffer_type());
        chunk_size_ = kMinChunkSize;
        omp_unset_lock(&lock_alloc_);
    }

  private:
    Pool(const pool_type &);
    const pool_type &operator =(const pool_type &);

    std::vector<pointer> heads_;
    std::vector<buffer_type> buffers_;
    std::deque<chunk_type> chunks_;
    uint32_t chunk_size_;
    omp_lock_t lock_alloc_;
    allocator_type alloc_;
};

namespace std {
template <typename T> inline void swap(Pool<T> &x, Pool<T> &y) {
    x.swap(y);
}
}

#endif

