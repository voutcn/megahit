/*
 *  io-utility.h
 *  This file is a part of MEGAHIT
 *  
 *  Copyright (C) 2014 The University of Hong Kong
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IO_UTILITY_H_
#define IO_UTILITY_H_

#include <assert.h>
#include <string>
#include <vector>
#include <zlib.h>
#include "definitions.h"
#include "kseq.h"
#include "mem_file_checker-inl.h"

KSEQ_INIT(gzFile, gzread)

struct ContigPackage {
    const static unsigned kMaxNumChars = (1 << 28); // tunable

    std::string seqs;
    std::vector<int> seq_lengths;
    std::vector<int> start_pos;
    std::vector<multi_t> multiplicity;
    int cur_pos;

    void ReadContigs(kseq_t *seq, char *dna_map) {
        clear();
        while (kseq_read(seq) >= 0) {
            if (seq->seq.l == 0) {
                continue;
            }

            start_pos.push_back(seqs.size());
            
            for (unsigned i = 0; i < seq->seq.l; ++i) {
                seqs.push_back(dna_map[uint8_t(seq->seq.s[i])]);
            }

            seq_lengths.push_back(seq->seq.l);
            if (seqs.length() >= kMaxNumChars) {
                break;
            }
        }
    }

    void ReadMultiplicity(gzFile multi_file) {
        multiplicity.resize(size());
        int num_bytes = gzread(multi_file, &multiplicity[0], sizeof(multi_t) * size());
        assert((unsigned)num_bytes == sizeof(multi_t) * size());
    }

    void clear() {
        seqs.clear();
        seq_lengths.clear();
        start_pos.clear();
    }

    size_t size() {
        return seq_lengths.size();
    }

    char CharAt(int i, int j) {
        return seqs[start_pos[i] + j];
    }

    std::string GetString(int i) {
        return seqs.substr(start_pos[i], seq_lengths[i]);
    }
};

struct ReadPackage {
    const static int kMaxNumReads = (1 << 22); // tunable
    const static int kMaxBytePerBatch = (1 << 28); // 256M, tunable
    int reads_per_batch;

    ReadPackage(int max_read_len) {
        init(max_read_len);
    }

    ReadPackage() {
        packed_reads = NULL;
        max_read_len = 0;
    }

    void init(int max_read_len) {
        this->max_read_len = max_read_len;
        num_of_reads = 0;
        offset_bits = 1;
        while ((1 << offset_bits) - 1 < max_read_len) {
            ++offset_bits;
        }
        words_per_read = (max_read_len * 2 + offset_bits + 31) / 32;
        printf("words per read: %d\n", words_per_read);
        reads_per_batch = kMaxBytePerBatch / (words_per_read * 4);
        reads_per_batch = reads_per_batch > kMaxNumReads ? kMaxNumReads : reads_per_batch;
        packed_reads = (edge_word_t*) MallocAndCheck(sizeof(edge_word_t) * words_per_read * reads_per_batch, __FILE__, __LINE__);
        assert(packed_reads != NULL);
    }

    ~ReadPackage() {
        if (packed_reads != NULL) {
            free(packed_reads);
        }
    }

    void ReadFastxReads(kseq_t *seq, char *dna_map) {
        clear();
        edge_word_t *cur_p = packed_reads;
        while (kseq_read(seq) >= 0) {
            if (int(seq->seq.l) > max_read_len) {
                fprintf(stderr, "Warning, a read longer than max_read_len: %d > %d\n", (int)(seq->seq.l), max_read_len);
                continue;
            }

            edge_word_t w = 0;
            int index = 0;
            int cur_word = 0;
            for (unsigned i = 0; i < seq->seq.l; ++i) {
                w = (w << 2) | dna_map[uint8_t(seq->seq.s[i])];
                ++index;
                if (index % 16 == 0) {
                    *(cur_p + cur_word) = w;
                    w = 0;
                    index = 0;
                    cur_word++;
                }
            }
            if (index != 0) {
                *(cur_p + cur_word) = w << (32 - index * 2);
                ++cur_word;
            }
            while (cur_word < words_per_read) {
                *(cur_p + cur_word) = 0;
                ++cur_word;
            }
            *(cur_p + words_per_read - 1) |= seq->seq.l;

            num_of_reads++;
            cur_p += words_per_read;
            if (num_of_reads >= reads_per_batch) {
                break;
            }
        }
    }

    void ReadFastxReadsAndReverse(kseq_t *seq, char *dna_map) {
        clear();
        edge_word_t *cur_p = packed_reads;
        while (kseq_read(seq) >= 0) {
            if ((int)(seq->seq.l) > max_read_len) {
                fprintf(stderr, "Warning, a read longer than max_read_len: %d > %d\n", (int)(seq->seq.l), max_read_len);
                continue;
            }

            edge_word_t w = 0;
            int index = 0;
            int cur_word = 0;
            for (int i = seq->seq.l - 1; i >= 0; --i) {
                w = (w << 2) | dna_map[uint8_t(seq->seq.s[i])];
                ++index;
                if (index % 16 == 0) {
                    *(cur_p + cur_word) = w;
                    w = 0;
                    index = 0;
                    cur_word++;
                }
            }
            if (index != 0) {
                *(cur_p + cur_word) = w << (32 - index * 2);
                ++cur_word;
            }
            while (cur_word < words_per_read) {
                *(cur_p + cur_word) = 0;
                ++cur_word;
            }
            *(cur_p + words_per_read - 1) |= seq->seq.l;

            num_of_reads++;
            cur_p += words_per_read;
            if (num_of_reads >= reads_per_batch) {
                break;
            }
        }
    }

    void ReadBinaryReads(gzFile read_file) { 
        int num_bytes = gzread(read_file, packed_reads, 
                                  sizeof(edge_word_t) * words_per_read * reads_per_batch);
        assert(num_bytes % (words_per_read * sizeof(edge_word_t)) == 0); 
        num_of_reads = num_bytes / (words_per_read * sizeof(edge_word_t));
    }

    void clear() {
        num_of_reads = 0;
    }

    int length(int i) {
        return packed_reads[(int64_t)(i + 1) * words_per_read - 1] & ((1 << offset_bits) - 1);
    }

    edge_word_t* GetReadPtr(int64_t read_id) {
        return packed_reads + read_id * words_per_read;
    }

    uint8_t CharAt(int i, int j) {
        return (packed_reads[(int64_t)i * words_per_read + j / 16] >> (15 - j % 16) * 2) & 3;
    }

    int max_read_len;
    int64_t num_of_reads;
    int words_per_read;
    int offset_bits;
    edge_word_t *packed_reads;
};

struct EdgeReader {
    std::vector<gzFile> edges_files;
    std::vector<edge_word_t*> cur_p;
    std::vector<int> num_edges_in_buffer;
    std::vector<int> edge_idx;
    int kmer_k, words_per_edge, words_per_true_edge;
    static const int kBufferSize = 4096;
    edge_word_t *buffer;

    void init(const std::string &prefix, int num_files) {
        edges_files.resize(num_files);
        cur_p.resize(num_files);
        num_edges_in_buffer.resize(num_files);
        edge_idx.resize(num_files);

        for (int i = 0; i < (int)edges_files.size(); ++i) {
            static char file_name[10240];
            sprintf(file_name, "%s.%d", prefix.c_str(), i);
            edges_files[i] = gzopen(file_name, "r");
            assert(edges_files[i] != NULL);
        }
        gzread(edges_files[0], &kmer_k, sizeof(edge_word_t));
        gzread(edges_files[0], &words_per_edge, sizeof(edge_word_t));

        words_per_true_edge = ((kmer_k + 1) * 2 + 31) / 32;

        buffer = (edge_word_t*) MallocAndCheck(sizeof(edge_word_t) * kBufferSize * words_per_edge * num_files, __FILE__, __LINE__);
        assert(buffer != NULL);

        for (int i = 0; i < (int)edges_files.size(); ++i) {
            if (!Refill_(i)) {
                --i;
            }
        }
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

        buffer = (edge_word_t*) MallocAndCheck(sizeof(edge_word_t) * kBufferSize * words_per_edge * num_files, __FILE__, __LINE__);
        assert(buffer != NULL);

        for (int i = 0; i < (int)edges_files.size(); ++i) {
            if (!Refill_(i)) {
                --i;
            }
        }
    }

    void destroy() {
        for (int i = 0; i < (int)edges_files.size(); ++i) {
            gzclose(edges_files[i]);
        }
        free(buffer);
    }

    bool Refill_(int i) {
        cur_p[i] = buffer + i * kBufferSize * words_per_edge;
        int num_bytes = gzread(edges_files[i], cur_p[i], sizeof(edge_word_t) * words_per_edge * kBufferSize);
        num_edges_in_buffer[i] = num_bytes / (sizeof(edge_word_t) * words_per_edge);
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

    bool NextEdge(edge_word_t *edge_p) {
        if (edges_files.size() == 0) { return false; }
        int max_idx = 0;
        for (int i = 1; i < (int)edges_files.size(); ++i) {
            if (compare_(i, max_idx)) {
                max_idx = i;
            }
        }

        memcpy(edge_p, cur_p[max_idx], sizeof(edge_word_t) * words_per_edge);
        cur_p[max_idx] += words_per_edge;
        ++edge_idx[max_idx];
        if (edge_idx[max_idx] == num_edges_in_buffer[max_idx]) {
            Refill_(max_idx);
        }

        return true;
    }

    bool NextEdgeUnsorted(edge_word_t *edge_p) {
        if (edges_files.size() == 0) { return false; }

        memcpy(edge_p, cur_p[0], sizeof(edge_word_t) * words_per_edge);
        cur_p[0] += words_per_edge;
        ++edge_idx[0];
        if (edge_idx[0] == num_edges_in_buffer[0]) {
            Refill_(0);
        }

        return true;
    }
};

#endif // IO_UTILITY_H_