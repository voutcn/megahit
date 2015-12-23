#ifndef SDBG_MULTI_IO__
#define SDBG_MULTI_IO__

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <algorithm>

#include "definitions.h"
#include "utils.h"
#include "mem_file_checker-inl.h"

struct SdbgPartitionRecord {
    int thread_id;
    long long starting_offset;
    long long num_items;
    long long num_tips;
    long long num_large_mul;
    long long num_last1;
    long long num_w[9];

    SdbgPartitionRecord(): thread_id(-1), starting_offset(0), num_items(0),
        num_tips(0), num_large_mul(0), num_last1(0) {
        memset(num_w, 0, sizeof(num_w));
    }
};

class SdbgWriter {
  private:
    std::string file_prefix_;
    int num_threads_;
    int num_buckets_;

    std::vector<FILE *> files_;
    std::vector<int> cur_bucket_;
    std::vector<int64_t> cur_thread_offset_;	// offset in BYTE
    std::vector<SdbgPartitionRecord> p_rec_;

    bool is_opened_;
    int kmer_size_;
    int words_per_tip_label_;

  public:

    SdbgWriter() {}
    ~SdbgWriter() {
        destroy();
    }

    void set_num_threads(int num_threads) {
        num_threads_ = num_threads;
    }
    void set_file_prefix(const std::string &file_prefix) {
        file_prefix_ = file_prefix;
    }
    void set_kmer_size(int k) {
        kmer_size_ = k;
        words_per_tip_label_ = DivCeiling(k * 2, 32);
    }
    void set_num_buckets(int num_buckets) {
        num_buckets_ = num_buckets;
    }

    void init_files() {
        files_.resize(num_threads_);
        cur_bucket_.resize(num_threads_, -1);
        cur_thread_offset_.resize(num_threads_, 0);
        p_rec_.resize(num_buckets_);

        for (int i = 0; i < num_threads_; ++i) {
            files_[i] = OpenFileAndCheck(FormatString("%s.sdbg.%d", file_prefix_.c_str(), i), "wb");
        }

        is_opened_ = true;
    }

    void write(int tid, int32_t bucket, int w, int last, int tip, multi_t multiplicity, uint32_t *packed_tip_label) {
        assert(tid < num_threads_);

        if (bucket != cur_bucket_[tid]) {
            cur_bucket_[tid] = bucket;
            assert(p_rec_[bucket].thread_id == -1);
            p_rec_[bucket].thread_id = tid;
            p_rec_[bucket].starting_offset = cur_thread_offset_[tid];
        }

        uint16_t packed_sdbg_item = w | (last << 4) | (tip << 5) | (std::min(multiplicity, (multi_t)kMulti2Sp) << 8);
        fwrite(&packed_sdbg_item, sizeof(uint16_t), 1, files_[tid]);
        ++p_rec_[bucket].num_items;
        ++p_rec_[bucket].num_w[w];
        p_rec_[bucket].num_last1 += last;
        cur_thread_offset_[tid] += sizeof(uint16_t);

        if (multiplicity > kMaxMulti2_t) {
            fwrite(&multiplicity, sizeof(multi_t), 1, files_[tid]);
            multiplicity = kMulti2Sp;
            ++p_rec_[bucket].num_large_mul;
            cur_thread_offset_[tid] += sizeof(multi_t);
        }

        if (tip) {
            fwrite(packed_tip_label, sizeof(uint32_t), words_per_tip_label_, files_[tid]);
            ++p_rec_[bucket].num_tips;
            cur_thread_offset_[tid] += sizeof(uint32_t) * words_per_tip_label_;
        }
    }

    int64_t num_edges() {
        int64_t total_edges = 0;

        for (int i = 0; i < num_buckets_; ++i) {
            total_edges += p_rec_[i].num_items;
        }

        return total_edges;
    }

    int64_t num_w(int w) {
        int64_t ret = 0;

        for (int i = 0; i < num_buckets_; ++i) {
            ret += p_rec_[i].num_w[w];
        }

        return ret;
    }

    int64_t num_last1() {
        int64_t ret = 0;

        for (int i = 0; i < num_buckets_; ++i) {
            ret += p_rec_[i].num_last1;
        }

        return ret;
    }

    int64_t num_tips() {
        int64_t ret = 0;

        for (int i = 0; i < num_buckets_; ++i) {
            ret += p_rec_[i].num_tips;
        }

        return ret;
    }

    void destroy() {
        if (is_opened_) {
            for (int i = 0; i < num_threads_; ++i) {
                fclose(files_[i]);
            }

            FILE *sdbg_info = OpenFileAndCheck((file_prefix_ + ".sdbg_info").c_str(), "w");
            fprintf(sdbg_info, "k %d\n", kmer_size_);
            fprintf(sdbg_info, "words_per_tip_label %d\n", words_per_tip_label_);
            fprintf(sdbg_info, "num_buckets %d\n", num_buckets_);
            fprintf(sdbg_info, "num_threads %d\n", num_threads_);

            int64_t total_edges = 0;
            int64_t total_tips = 0;
            int64_t total_large_mul = 0;

            for (int i = 0; i < num_buckets_; ++i) {
                total_edges += p_rec_[i].num_items;
                total_tips += p_rec_[i].num_tips;
                total_large_mul += p_rec_[i].num_large_mul;
            }

            fprintf(sdbg_info, "total_size %lld\n", (long long)total_edges);
            fprintf(sdbg_info, "num_tips %lld\n", (long long)total_tips);
            fprintf(sdbg_info, "large_multi %lld\n", (long long)total_large_mul);

            for (int i = 0; i < num_buckets_; ++i) {
                fprintf(sdbg_info, "%d %d %lld %lld %lld %lld\n", i,
                        p_rec_[i].thread_id,
                        (long long)p_rec_[i].starting_offset,
                        (long long)p_rec_[i].num_items,
                        (long long)p_rec_[i].num_tips,
                        (long long)p_rec_[i].num_large_mul);
            }

            fclose(sdbg_info);

            files_.clear();
            cur_bucket_.clear();
            cur_thread_offset_.clear();	// offset in BYTE
            p_rec_.clear();

            is_opened_ = false;
        }
    }
};

class SdbgReader {
  private:
    std::string file_prefix_;
    int kmer_size_;
    int words_per_tip_label_;
    long long num_items_;
    long long num_tips_;
    long long num_large_mul_;
    long long f_[6];
    int num_files_;
    int num_buckets_;

    int pre_lkt_len_;
    int pre_lkt_size_;

    int cur_bucket_;
    long long cur_bucket_cnt_;
    long long cur_vol_;
    void *cur_bucket_ptr_;

    std::vector<SdbgPartitionRecord> p_rec_;
    std::vector<long long> pre_lkt_;
    std::vector<long long> file_sizes_;
    std::vector<int> fds_;
    void *mmap_;
    long page_size_;
    int64_t mmap_size_;

    bool is_opened_;

  public:
    SdbgReader(): is_opened_(false) {}
    ~SdbgReader() {
        destroy();
    }

    void set_file_prefix(const std::string &prefix) {
        file_prefix_ = prefix;
    }

    void read_info() {
        FILE *sdbg_info = OpenFileAndCheck((file_prefix_ + ".sdbg_info").c_str(), "r");
        assert(fscanf(sdbg_info, "k %d\n", &kmer_size_) == 1);
        assert(fscanf(sdbg_info, "words_per_tip_label %d\n", &words_per_tip_label_) == 1);
        assert(fscanf(sdbg_info, "num_buckets %d\n", &num_buckets_) == 1);
        assert(fscanf(sdbg_info, "num_threads %d\n", &num_files_) == 1);
        assert(fscanf(sdbg_info, "total_size %lld\n", &num_items_) == 1);
        assert(fscanf(sdbg_info, "num_tips %lld\n", &num_tips_) == 1);
        assert(fscanf(sdbg_info, "large_multi %lld\n", &num_large_mul_) == 1);

        p_rec_.resize(num_buckets_);
        file_sizes_.resize(num_files_, 0);
        long long acc = 0;
        f_[0] = -1;
        f_[1] = 0;

        for (int i = 0; i < num_buckets_; ++i) {
            int dummy;
            assert(fscanf(sdbg_info, "%d %d %lld %lld %lld %lld\n", &dummy,
                          &p_rec_[i].thread_id,
                          &p_rec_[i].starting_offset,
                          &p_rec_[i].num_items,
                          &p_rec_[i].num_tips,
                          &p_rec_[i].num_large_mul) == 6);
            file_sizes_[p_rec_[i].thread_id] += p_rec_[i].num_items * sizeof(uint16_t) + p_rec_[i].num_tips * sizeof(uint32_t) * words_per_tip_label_ + p_rec_[i].num_large_mul * sizeof(multi_t);
            acc += p_rec_[i].num_items;
            f_[i / (num_buckets_ / 4) + 2] = acc;
        }

        if (fscanf(sdbg_info, "%d", &pre_lkt_len_) == EOF) {
            pre_lkt_len_ = 0;
            pre_lkt_size_ = num_buckets_;
            while ((1 << (pre_lkt_len_ * 2)) < num_buckets_) {
                ++pre_lkt_len_;
            }
            pre_lkt_.resize(pre_lkt_size_ * 2);

            for (long long acc = 0, i = 0; i < pre_lkt_size_; ++i) {
                pre_lkt_[i * 2] = acc;
                acc += p_rec_[i].num_items;
                pre_lkt_[i * 2 + 1] = acc - 1;
            }
        } else {
            pre_lkt_size_ = 1 << (pre_lkt_len_ * 2);
            pre_lkt_.resize(pre_lkt_size_ * 2);
            for (int i = 0; i < pre_lkt_size_; ++i) {
                assert(fscanf(sdbg_info, "%lld%lld", &pre_lkt_[i*2], &pre_lkt_[i*2+1]) == 2);
            }
        }

        fclose(sdbg_info);
    }

    int kmer_size() const {
        return kmer_size_;
    }
    int words_per_tip_label() const {
        return words_per_tip_label_;
    }
    int prefix_lkt_len() const {
        return pre_lkt_len_;
    }
    int prefix_lkt_size() const {
        return pre_lkt_size_;
    }
    long long num_items() const {
        return num_items_;
    }
    long long num_tips() const {
        return num_tips_;
    }
    long long num_large_mul() const {
        return num_large_mul_;
    }
    const long long *f() const {
        return f_;
    }
    long long bucket_size(int i) const {
        return p_rec_[i].num_items;
    }
    long long prefix_lkt(int i) const {
        return pre_lkt_[i];
    }

    void init_files() {
        assert(!is_opened_);
        fds_.resize(num_files_);
        page_size_ = sysconf(_SC_PAGESIZE);

        for (int i = 0; i < num_files_; ++i) {
            fds_[i] = open(FormatString("%s.sdbg.%d", file_prefix_.c_str(), i), O_RDONLY);
            assert(fds_[i] != -1);
        }

        cur_bucket_ = -1;
        cur_bucket_cnt_ = 0;
        cur_vol_ = 0;
        mmap_ = NULL;
        mmap_size_ = 0;
        is_opened_ = true;
    }

    bool NextItem(uint16_t &item) {
        if (cur_bucket_ >= num_buckets_) {
            return false;
        }

        while (cur_bucket_cnt_ >= cur_vol_) {
            cur_bucket_cnt_ = 0;
            ++cur_bucket_;

            while (cur_bucket_ < num_buckets_ && p_rec_[cur_bucket_].thread_id == -1) {
                ++cur_bucket_;
            }

            if (cur_bucket_ >= num_buckets_) {
                return false;
            }

            if (mmap_) {
                munmap(mmap_, mmap_size_);
                mmap_ = NULL;
            }

            int64_t offset = p_rec_[cur_bucket_].starting_offset / page_size_ * page_size_;
            mmap_size_ = p_rec_[cur_bucket_].num_items * sizeof(uint16_t) +
                         p_rec_[cur_bucket_].num_tips * sizeof(uint32_t) * words_per_tip_label_ +
                         p_rec_[cur_bucket_].num_large_mul * sizeof(multi_t);
            mmap_size_ += p_rec_[cur_bucket_].starting_offset - offset;

            mmap_ = mmap(NULL, mmap_size_, PROT_READ, MAP_PRIVATE, fds_[p_rec_[cur_bucket_].thread_id], offset);
            assert(mmap_ != NULL);
            madvise(mmap_, mmap_size_, MADV_SEQUENTIAL);

            cur_bucket_ptr_ = (char*)mmap_ + p_rec_[cur_bucket_].starting_offset - offset;
            cur_vol_ = p_rec_[cur_bucket_].num_items;
        }

        ++cur_bucket_cnt_;
        item = *((uint16_t *)cur_bucket_ptr_);
        cur_bucket_ptr_ = (void *)((char *)cur_bucket_ptr_ + sizeof(uint16_t));
        return true;
    }

    multi_t NextLargeMul() {
        multi_t large_mul = *((multi_t *)cur_bucket_ptr_);
        cur_bucket_ptr_ = (void *)((char *)cur_bucket_ptr_ + sizeof(multi_t));
        return large_mul;
    }

    void NextTipLabel(uint32_t *tip_label) {
        memcpy(tip_label, cur_bucket_ptr_, sizeof(uint32_t) * words_per_tip_label_);
        cur_bucket_ptr_ = (void *)((char *)cur_bucket_ptr_ + sizeof(uint32_t) * words_per_tip_label_);
    }

    void destroy() {
        if (is_opened_) {
            if (mmap_) {
                munmap(mmap_, mmap_size_);
                mmap_ = NULL;
            }
            for (int i = 0; i < num_files_; ++i) {
                close(fds_[i]);
            }

            file_sizes_.clear();
            fds_.clear();

            is_opened_ = false;
        }
    }
};

#endif