/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong
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

#ifndef SEQUENCE_MANAGER_H__
#define SEQUENCE_MANAGER_H__

#include <string>
#include <vector>

#include <zlib.h>

#include "definitions.h"
#include "kseq.h"
#include "sequence_package.h"
#include "edge_io.h"

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

/**
 * @brief   manage the reading of a set of fastx reads/binary reads/contigs/edges
 * @details providing functions like reading contigs/single-end reads/paired reads/interleaved reads from fastx/binary format
 * 			and writing sequences into binary formats
 */
struct SequenceManager {
    enum FileType {
        kFastxReads,
        kBinaryReads,
        kMegahitContigs,
        kMegahitEdges,
        kSortedEdges,
    } f_type;

    enum ReadLibType {
        kPaired,
        kInterleaved,
        kSingle,
        kOthers,
    } r_type;

    SequencePackage *package_;
    std::vector<multi_t> *multi_; // edges or contigs' multiplicity
    std::vector<float> *float_multi_;
    std::vector<gzFile> files_; // for reading sorted edges
    std::vector<kseq_t *> kseq_readers_; // for paired reading
    EdgeReader edge_reader_;
    bool edge_reader_inited_;

    std::vector<uint32_t> buf_;		// for reading binary reads
    int min_len_;
    int k_from_, k_to_;

    SequenceManager(SequencePackage *package = NULL): package_(package) {
        multi_ = NULL;
        float_multi_ = NULL;
        files_.clear();
        kseq_readers_.clear();
        edge_reader_.destroy();
        edge_reader_inited_ = false;
        min_len_ = 0;
        k_from_ = 1;
        k_to_ = 1;
    };

    ~SequenceManager() {
        clear();
    }

    void clear() {
        for (auto iter = kseq_readers_.begin(); iter != kseq_readers_.end(); ++iter) {
            kseq_destroy(*iter);
        }

        for (auto iter = files_.begin(); iter != files_.end(); ++iter) {
            gzclose(*iter);
        }

        files_.clear();
        kseq_readers_.clear();

        if (edge_reader_inited_) {
            edge_reader_.destroy();
            edge_reader_inited_ = false;
        }
    }

    void set_file_type(FileType f_type) {
        this->f_type = f_type;
    }
    void set_readlib_type(ReadLibType r_type) {
        this->r_type = r_type;
    }
    void set_package(SequencePackage *package) {
        package_ = package;
    }
    void set_multiplicity_vector(std::vector<multi_t> *v) {
        multi_ = v;
    }
    void set_float_multiplicity_vector(std::vector<float> *v) {
        float_multi_ = v;
    }

    void set_file(const std::string &file_name);
    void set_pe_files(const std::string &file_name1, const std::string &file_name2);
    void set_edge_files(const std::string &file_prefix);
    void set_kmer_size(int kmer_from, int kmer_to) {
        k_from_ = kmer_from;    // only apply to reading megahit contigs
        k_to_ = kmer_to;
    }
    void set_min_len(int min_len) {
        min_len_ = min_len;    // only apply to megahit contigs
    }

    int64_t ReadShortReads(int64_t max_num, int64_t max_num_bases, bool append, bool reverse, bool trimN = false, std::string file_name = "");
    int64_t ReadEdges(int64_t max_num, bool append);
    int64_t ReadEdgesWithFixedLen(int64_t max_num, bool append);
    int64_t ReadShortedEdges(int64_t max_num, int kmer_size, bool append);
    int64_t ReadMegahitContigs(int64_t max_num, int64_t max_num_bases, bool append, bool reverse,
                               int discard_flag, bool extend_loop, bool calc_depth);
    void WriteBinarySequences(FILE *file, bool reverse, int64_t from = 0, int64_t to = -1);
};

#endif