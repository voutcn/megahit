/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
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

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#ifndef LOCAL_ASSEMBLER_H__
#define LOCAL_ASSEMBLER_H__

#include <stdint.h>
#include <vector>
#include <string>
#include <deque>

#include <omp.h>

#include "hash_table.h"
#include "sequence_package.h"
#include "kmer_plus.h"
#include "lib_info.h"

struct LocalAssembler {
    static const int kMaxLocalRange = 650;
    static const int kMaxNumLocks = 1 << 20;

    typedef std::pair<double, double> tlen_t;
    typedef KmerPlus<2, uint32_t, uint64_t> kmer_plus_t;
    typedef Kmer<2, uint32_t> kmer_t;
    typedef HashTable<kmer_plus_t, kmer_t> mapper_t;

    struct MappingRecord {
        uint32_t contig_id;
        int32_t contig_from : 28;
        int32_t contig_to : 28;
        int16_t query_from : 15;
        int16_t query_to : 15;
        bool strand : 1;
        int mismatch: 9;
    };

    int min_contig_len_;	// only align reads to these contigs
    int seed_kmer_;			// kmer size for seeding
    int similarity_;		// similarity threshold for alignment
    int sparsity_;			// sparsity of hash mapper
    int min_mapped_len_;
    std::string local_filename_;
    int local_kmin_, local_kmax_, local_step_;
    int num_threads_;

    // auto calculated after calling EstimateInsertSize() and MapToContigs()
    int local_range_;
    int max_read_len_;

    SequencePackage *contigs_;
    SequencePackage *reads_;
    mapper_t mapper_;
    std::vector<lib_info_t> lib_info_;
    std::vector<tlen_t> insert_sizes_;
    std::vector<omp_lock_t> locks_;
    volatile int output_lock_;

    std::vector<std::deque<uint64_t> > mapped_f_, mapped_r_;

    LocalAssembler(int min_contig_len, int seed_kmer, int sparsity)
        : min_contig_len_(min_contig_len), seed_kmer_(seed_kmer), sparsity_(sparsity),
          contigs_(NULL), reads_(NULL) {
        similarity_ = 0.95;
        min_mapped_len_ = 100;
        num_threads_ = omp_get_max_threads();

        local_kmin_ = 21;
        local_kmax_ = 41;
        local_step_ = 6;
    }

    ~LocalAssembler() {
        delete contigs_;
        delete reads_;
    }

    void set_kmer(int kmin, int kmax, int step) {
        local_kmin_ = kmin;
        local_kmax_ = kmax;
        local_step_ = step;
    }

    void set_mapping_threshold(double similarity, int mapping_len) {
        similarity_ = similarity;
        min_mapped_len_ = mapping_len;
    }

    void set_num_threads(int num_threads) {
        num_threads_ = num_threads;
    }

    void set_local_file(const std::string &local_filename) {
        local_filename_ = local_filename;
    }

    void ReadContigs(const std::string &file_name);
    void BuildHashMapper(bool show_stat = true);
    void AddReadLib(const std::string &file_prefix);
    void EstimateInsertSize(bool show_stat = true);
    void MapToContigs();
    void LocalAssemble();

    void AddToHashMapper_(mapper_t &mapper, unsigned contig_id, int sparsity);
    int Match_(size_t read_id, int query_from, int query_to, size_t contig_id, int ref_from, int ref_to, bool strand);
    int LocalRange_(int lib_id);
    int AddToMappingDeque_(size_t read_id, const MappingRecord &rec, int local_range);
    int AddMateToMappingDeque_(size_t read_id, size_t mate_id, const MappingRecord &rec1, const MappingRecord &rec2, bool mapped2, int local_range);
    bool MapToHashMapper_(const mapper_t &mapper, size_t read_id, MappingRecord &rec);
    static void *LocalAssembleThread_(void *task);
};

#endif