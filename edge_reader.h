#ifndef EDGE_READER_H__
#define EDGE_READER_H__

#include <zlib.h>
#include <assert.h>
#include <vector>

#include "utils.h"
#include "mem_file_checker-inl.h"
#include "definitions.h"

struct EdgeReader2 {
    std::vector<gzFile> edges_files;
    std::vector<uint32_t*> cur_p;
    std::vector<int> num_edges_in_buffer;
    std::vector<int> edge_idx;
    int kmer_k, words_per_edge, words_per_true_edge;
    static const int kBufferSize = 4096;
    uint32_t *buffer;
    int min_idx_;

    void init(const std::string &prefix, int num_files) {
        edges_files.resize(num_files);
        cur_p.resize(num_files);
        num_edges_in_buffer.resize(num_files);
        edge_idx.resize(num_files);

        for (int i = 0; i < (int)edges_files.size(); ++i) {
            edges_files[i] = gzopen(FormatString("%s.%d", prefix.c_str(), i), "r");
            assert(edges_files[i] != NULL);
        }
        gzread(edges_files[0], &kmer_k, sizeof(uint32_t));
        gzread(edges_files[0], &words_per_edge, sizeof(uint32_t));

        words_per_true_edge = ((kmer_k + 1) * 2 + 31) / 32;

        buffer = (uint32_t*) MallocAndCheck(sizeof(uint32_t) * kBufferSize * words_per_edge * num_files, __FILE__, __LINE__);
        assert(buffer != NULL);

        for (int i = 0; i < (int)edges_files.size(); ++i) {
            if (!Refill_(i)) {
                --i;
            }
        }
        min_idx_ = -1;
    }

    void InitUnsorted(const std::string &prefix, int num_files, int kmer_k, int words_per_edge) {
        edges_files.resize(num_files);
        cur_p.resize(num_files);
        num_edges_in_buffer.resize(num_files);
        edge_idx.resize(num_files);

        this->words_per_edge = words_per_edge;
        this->kmer_k = kmer_k;

        for (int i = 0; i < (int)edges_files.size(); ++i) {
            static char file_name[10240];
            sprintf(file_name, "%s.%d", prefix.c_str(), i);
            edges_files[i] = gzopen(file_name, "r");
            assert(edges_files[i] != NULL);
        }

        words_per_true_edge = ((kmer_k + 1) * 2 + 31) / 32;

        buffer = (uint32_t*) MallocAndCheck(sizeof(uint32_t) * kBufferSize * words_per_edge * num_files, __FILE__, __LINE__);
        assert(buffer != NULL);

        for (int i = 0; i < (int)edges_files.size(); ++i) {
            if (!Refill_(i)) {
                --i;
            }
        }
        min_idx_ = -1;
    }

    void destroy() {
        for (int i = 0; i < (int)edges_files.size(); ++i) {
            gzclose(edges_files[i]);
        }
        free(buffer);
    }

    bool Refill_(int i) {
        cur_p[i] = buffer + i * kBufferSize * words_per_edge;
        int num_bytes = gzread(edges_files[i], cur_p[i], sizeof(uint32_t) * words_per_edge * kBufferSize);
        num_edges_in_buffer[i] = num_bytes / (sizeof(uint32_t) * words_per_edge);
        edge_idx[i] = 0;
        if (num_edges_in_buffer[i] == 0) {
            gzclose(edges_files[i]);
            edges_files.erase(edges_files.begin() + i);
            cur_p.erase(cur_p.begin() + i);
            num_edges_in_buffer.erase(num_edges_in_buffer.begin() + i);
            edge_idx.erase(edge_idx.begin() + i);

            return false;
        }

        return true;
    }

    bool compare_(int i, int j) {
        for (int k = 0; k < words_per_true_edge; ++k) {
            if (cur_p[i][k] < cur_p[j][k]) return true;
            if (cur_p[i][k] > cur_p[j][k]) return false;
        }
        assert(false);
    }

    uint32_t* NextSortedEdge() {
        if (edges_files.size() == 0) {
            return NULL;
        }

        if(min_idx_ >= 0) {
            cur_p[min_idx_] += words_per_edge;
            ++edge_idx[min_idx_];
            if (edge_idx[min_idx_] == num_edges_in_buffer[min_idx_]) {
                Refill_(min_idx_);
            }
        }

        if (edges_files.size() == 0) {
            return NULL;
        }

        for (int i = 1; i < (int)edges_files.size(); ++i) {
            if (compare_(i, min_idx_)) {
                min_idx_ = i;
            }
        }

        return cur_p[min_idx_];
    }

    uint32_t* NextUnsortedEdge() {
        if (edges_files.size() == 0) {
            return NULL;
        }

        if (min_idx_ != -1) {
            cur_p[0] += words_per_edge;
            ++edge_idx[0];
            if (edge_idx[0] == num_edges_in_buffer[0]) {
                Refill_(0);
            }
        } else {
            min_idx_ = 0;
        }

        if (edges_files.size() == 0) {
            return NULL;
        }

        return cur_p[0];
    }
};

#endif