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

#include "local_assembler.h"

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include <zlib.h>
#include <omp.h>
#include "lib_idba/sequence.h"
#include "lib_idba/hash_graph.h"
#include "lib_idba/contig_graph.h"

#include "utils.h"
#include "mem_file_checker-inl.h"
#include "kseq.h"
#include "histgram.h"
#include "bit_operation.h"
#include "sequence_manager.h"
#include "read_lib_functions-inl.h"

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

void LocalAssembler::ReadContigs(const std::string &contig_file_name) {
    SequenceManager seq_manager;
    seq_manager.set_file_type(SequenceManager::kMegahitContigs);
    seq_manager.set_file(contig_file_name);
    seq_manager.set_min_len(min_contig_len_);

    if (contigs_ == NULL) {
        contigs_ = new SequencePackage();
    }

    seq_manager.set_package(contigs_);
    bool contig_reverse = false;
    bool append_to_package = false;
    int discard_flag = contig_flag::kLoop;
    bool extend_loop = false;
    bool calc_depth = false;

    seq_manager.ReadMegahitContigs(1LL << 60, 1LL << 60, append_to_package, contig_reverse, discard_flag, extend_loop, calc_depth);
}

void LocalAssembler::BuildHashMapper(bool show_stat) {
    size_t sz = contigs_->size();
    size_t estimate_num_kmer = 0;

    #pragma omp parallel for reduction(+: estimate_num_kmer)

    for (size_t i = 0; i < sz; ++i) {
        estimate_num_kmer += (contigs_->length(i) - seed_kmer_ + sparsity_) / sparsity_;
    }

    mapper_.reserve(estimate_num_kmer);

    #pragma omp parallel for

    for (size_t i = 0; i < sz; ++i) {
        AddToHashMapper_(mapper_, i, sparsity_);
    }

    if (show_stat) {
        xlog("Number of contigs: %lu, Mapper size: %lu\n", contigs_->size(), mapper_.size());
    }
}

void LocalAssembler::AddReadLib(const std::string &file_prefix) {
    if (reads_ == NULL) {
        reads_ = new SequencePackage();
    }

    ReadBinaryLibs(file_prefix, *reads_, lib_info_);
    insert_sizes_.resize(lib_info_.size(), tlen_t(-1, -1));
}

inline uint64_t EncodeContigOffset(unsigned contig_id, unsigned contig_offset, bool strand) {
    return (uint64_t(contig_id) << 32) | (contig_offset << 1) | strand;
}

inline void DecodeContigOffset(uint64_t code, uint32_t &contig_id, uint32_t &contig_offset, bool &strand) {
    contig_id = code >> 32;
    contig_offset = (code & 0xFFFFFFFFULL) >> 1;
    strand = code & 1ULL;
}

void LocalAssembler::AddToHashMapper_(mapper_t &mapper, unsigned contig_id, int sparcity) {
    kmer_plus_t kp;
    kmer_t key;
    kp.ann = ~0ULL;	// special marker

    for (int i = 0, len = contigs_->length(contig_id); i + seed_kmer_ <= len; i += sparcity) {
        uint64_t full_offset = contigs_->get_start_index(contig_id) + i;
        key.init(&contigs_->packed_seq[full_offset / SequencePackage::kCharsPerWord], full_offset % SequencePackage::kCharsPerWord, seed_kmer_);
        kp.kmer = key.unique_format(seed_kmer_);
        kmer_plus_t &kp_in_mapper = mapper.find_or_insert_with_lock(kp);

        if (kp_in_mapper.ann != ~0ULL) {
            kp_in_mapper.ann |= 1ULL << 63; // mark the highest bit as unused
        }
        else {
            kp_in_mapper.ann = EncodeContigOffset(contig_id, i, key != kp_in_mapper.kmer);
        }

        mapper.unlock(kp);
    }
}

inline uint32_t GetWord(uint32_t *first_word, uint32_t first_shift, int from, int len, bool strand) {
    int from_word_idx = (first_shift + from) / 16;
    int from_word_shift = (first_shift + from) % 16;
    uint32_t ret = *(first_word + from_word_idx) << from_word_shift * 2;
    assert(len <= 16);

    if (16 - from_word_shift < len) {
        ret |= *(first_word + from_word_idx + 1) >> (16 - from_word_shift) * 2;
    }

    if (len < 16) {
        ret >>= (16 - len) * 2;
        ret <<= (16 - len) * 2;
    }

    if (strand == 1) {
        bit_operation::ReverseComplement(ret);
        ret <<= (16 - len) * 2;
    }

    return ret;
}

inline int Mismatch(uint32_t x, uint32_t y) {
    x ^= y;
    x |= x >> 1;
    x &= 0x55555555U;
    return __builtin_popcount(x);
}

int LocalAssembler::Match_(size_t read_id, int query_from, int query_to,
                           size_t contig_id, int ref_from, int ref_to, bool strand) {
    uint32_t *query_first_word = &reads_->packed_seq[reads_->get_start_index(read_id) / 16];
    int query_shift = reads_->get_start_index(read_id) % 16;
    uint32_t *ref_first_word = &contigs_->packed_seq[contigs_->get_start_index(contig_id) / 16];
    int ref_shift = contigs_->get_start_index(contig_id) % 16;

    int match_len = query_to - query_from + 1;
    int threshold = similarity_ * match_len + 0.5;

    for (int i = query_from; i <= query_to; i += 16) {
        int len = std::min(16, query_to - i + 1);
        uint32_t qw = GetWord(query_first_word, query_shift, i, len, 0);
        int ref_i = strand == 0 ? ref_from + i - query_from : ref_to - (i + len - 1 - query_from);
        uint32_t rw = GetWord(ref_first_word, ref_shift, ref_i, len, strand);

        match_len -= Mismatch(qw, rw);

        if (match_len < threshold) {
            return 0;
        }
    }

    return match_len;
}

bool LocalAssembler::MapToHashMapper_(const mapper_t &mapper, size_t read_id, MappingRecord &rec) {
    int len = reads_->length(read_id);

    if (len < seed_kmer_ || len < 50) return false; // too short reads not reliable

    int tested = 0;
    MappingRecord tested_rec[3];

    uint32_t *packed_seq = &reads_->packed_seq[0];
    kmer_t kmer_f(packed_seq + reads_->get_start_index(read_id) / SequencePackage::kCharsPerWord,
                  reads_->get_start_index(read_id) % SequencePackage::kCharsPerWord, seed_kmer_);
    kmer_t kmer_r = kmer_f;
    kmer_r.ReverseComplement(seed_kmer_);
    int num_mapped = 0;
    uint32_t contig_id, contig_offset;
    bool contig_strand;

    for (int i = seed_kmer_ - 1; i < len; ++i) {
        if (i >= seed_kmer_) {
            uint8_t ch = reads_->get_base(read_id, i);
            kmer_f.ShiftAppend(ch, seed_kmer_);
            kmer_r.ShiftPreappend(3 - ch, seed_kmer_);
        }

        bool query_strand = kmer_f.cmp(kmer_r, seed_kmer_) <= 0 ? 0 : 1;

        auto iter = mapper.find(query_strand == 0 ? kmer_f : kmer_r);

        if (iter == mapper.end() || (iter->ann >> 63) != 0) {
            continue;
        }

        DecodeContigOffset(iter->ann, contig_id, contig_offset, contig_strand);
        assert(contig_id < contigs_->size());
        assert(contig_offset < contigs_->length(contig_id));

        bool mapping_strand = contig_strand ^ query_strand;
        int contig_from = mapping_strand == 0 ? contig_offset - (i - seed_kmer_ + 1) : contig_offset - (len - 1 - i);
        int contig_to = mapping_strand == 0 ? contig_offset + seed_kmer_ - 1 + len - 1 - i : contig_offset + i;
        contig_from = std::max(contig_from, 0);
        contig_to = std::min((int)contigs_->length(contig_id) - 1, contig_to);

        if (contig_to - contig_from + 1 < len && contig_to - contig_from + 1 < min_mapped_len_) {
            continue;    // clipped alignment is considered iff its length >= min_mapped_len_
        }

        int query_from = mapping_strand == 0 ? i - (seed_kmer_ - 1) - (contig_offset - contig_from) : i - (contig_to - contig_offset);
        int query_to = mapping_strand == 0 ? i - (seed_kmer_ - 1) + (contig_to - contig_offset) : i + (contig_offset - contig_from);

        bool has_tested = false;

        for (int j = 0; j < tested; ++j) {
            if (contig_id == tested_rec[j].contig_id &&
                    contig_from == tested_rec[j].contig_from &&
                    contig_to == tested_rec[j].contig_to &&
                    query_from == tested_rec[j].query_from &&
                    query_to == tested_rec[j].query_to &&
                    mapping_strand == tested_rec[j].strand) {
                has_tested = true;
                break;
            }
        }

        if (has_tested) {
            continue;
        }
        else {
            if (tested >= 3) {
                tested--;
            }

            tested_rec[tested].contig_id = contig_id;
            tested_rec[tested].query_from = query_from;
            tested_rec[tested].query_to = query_to;
            tested_rec[tested].contig_from = contig_from;
            tested_rec[tested].contig_to = contig_to;
            tested_rec[tested].strand = mapping_strand;
            ++tested;
        }

        int match_bases = Match_(read_id, query_from, query_to, contig_id, contig_from, contig_to, mapping_strand);

        if (match_bases > 0) {
            if (num_mapped > 0) {
                return false;
            }
            else {
                rec = tested_rec[tested - 1];
                rec.mismatch = query_to - query_from + 1 - match_bases;
                num_mapped = 1;
            }
        }
    }

    return num_mapped == 1;
}

void LocalAssembler::EstimateInsertSize(bool show_stat) {
    for (unsigned lib_id = 0; lib_id < lib_info_.size(); ++lib_id) {
        if (!lib_info_[lib_id].is_pe) {
            continue;
        }

        MappingRecord rec1, rec2;
        Histgram<int> insert_hist;
        int64_t start_read_id = lib_info_[lib_id].from;
        int64_t end_read_id = start_read_id;

        while (insert_hist.size() < (1 << 18) && end_read_id <= lib_info_[lib_id].to) {
            start_read_id = end_read_id;
            end_read_id = std::min(lib_info_[lib_id].to + 1, start_read_id + (2 << 18));

            #pragma omp parallel for private(rec1, rec2)

            for (int64_t i = start_read_id; i < end_read_id; i += 2) {
                if (MapToHashMapper_(mapper_, i, rec1) &&
                        MapToHashMapper_(mapper_, i + 1, rec2)) {
                    if (rec1.contig_id == rec2.contig_id && rec1.strand != rec2.strand) {
                        int insert_size = -1;

                        if (rec1.strand == 0) {
                            insert_size = rec2.contig_to + reads_->length(i + 1) - rec2.query_to - (rec1.contig_from - rec1.query_from);
                        }
                        else {
                            insert_size = rec1.contig_to + reads_->length(i) - rec1.query_to - (rec2.contig_from - rec2.query_from);
                        }

                        if (insert_size >= (int)reads_->length(i) &&
                                insert_size >= (int)reads_->length(i + 1)) {
                            insert_hist.insert(insert_size);
                        }
                    }
                }
            }
        }

        insert_hist.Trim(0.01);
        insert_sizes_[lib_id] = tlen_t(insert_hist.mean(), insert_hist.sd());

        if (show_stat) {
            xlog("Lib %d, insert size: %.2lf sd: %.2lf\n", lib_id, insert_hist.mean(), insert_hist.sd());
        }
    }
}

int LocalAssembler::LocalRange_(int lib_id) {
    int local_range = lib_info_[lib_id].max_read_len - 1;

    if (insert_sizes_[lib_id].first >= lib_info_[lib_id].max_read_len) {
        local_range = std::min(2 * insert_sizes_[lib_id].first,
                               insert_sizes_[lib_id].first + 3 * insert_sizes_[lib_id].second);
    }

    if (local_range > kMaxLocalRange) {
        local_range = kMaxLocalRange;
    }

    return local_range;
}

inline uint64_t PackMappingResult(uint64_t contig_offset, uint64_t is_mate, uint64_t mismatch,
                                  uint64_t strand, uint64_t read_id) {
    assert(contig_offset < (1ULL << 14));
    return (contig_offset << 50) | (is_mate << 49) | (std::min(uint64_t(15), mismatch) << 45) |
           (strand << 44) | read_id;
}

int LocalAssembler::AddToMappingDeque_(size_t read_id, const MappingRecord &rec, int local_range) {
    // xlog("%lu %Lu\n", read_id, rec.contig_id);
    assert(read_id < reads_->size());
    assert(rec.contig_id < contigs_->size());

    int contig_len = contigs_->length(rec.contig_id);
    int read_len = reads_->length(read_id);
    int ret = 0;

    if (rec.contig_to < local_range && rec.query_from != 0 && rec.query_to == read_len - 1) {
        uint64_t res = PackMappingResult(rec.contig_to, 0, rec.mismatch, rec.strand, read_id);
        omp_set_lock(&locks_[rec.contig_id % kMaxNumLocks]);
        mapped_f_[rec.contig_id].push_back(res);
        omp_unset_lock(&locks_[rec.contig_id % kMaxNumLocks]);
        ret++;
    }
    else if (rec.contig_from + local_range >= contig_len && rec.query_to < read_len - 1 && rec.query_from == 0) {
        uint64_t res = PackMappingResult(contig_len - 1 - rec.contig_from, 0, rec.mismatch, rec.strand, read_id);
        omp_set_lock(&locks_[rec.contig_id % kMaxNumLocks]);
        mapped_r_[rec.contig_id].push_back(res);
        omp_unset_lock(&locks_[rec.contig_id % kMaxNumLocks]);
        ret++;
    }

    return ret;
}

int LocalAssembler::AddMateToMappingDeque_(size_t read_id, size_t mate_id, const MappingRecord &rec1, const MappingRecord &rec2, bool mapped2, int local_range) {
    assert(read_id < reads_->size());
    assert(mate_id < reads_->size());
    assert(rec1.contig_id < contigs_->size());
    assert(!mapped2 || rec2.contig_id < contigs_->size());

    if (mapped2 && rec2.contig_id == rec1.contig_id)
        return 0;

    int contig_len = contigs_->length(rec1.contig_id);
    int ret = 0;

    if (rec1.contig_to < local_range && rec1.strand == 1) {
        uint64_t res = PackMappingResult(rec1.contig_to, 1, rec1.mismatch, rec1.strand, mate_id);
        omp_set_lock(&locks_[rec1.contig_id % kMaxNumLocks]);
        mapped_f_[rec1.contig_id].push_back(res);
        omp_unset_lock(&locks_[rec1.contig_id % kMaxNumLocks]);
        ret++;
    }
    else if (rec1.contig_from + local_range >= contig_len && rec1.strand == 0) {
        uint64_t res = PackMappingResult(contig_len - 1 - rec1.contig_from, 1, rec1.mismatch, rec1.strand, mate_id);
        omp_set_lock(&locks_[rec1.contig_id % kMaxNumLocks]);
        mapped_r_[rec1.contig_id].push_back(res);
        omp_unset_lock(&locks_[rec1.contig_id % kMaxNumLocks]);
        ret++;
    }

    return ret;
}

void LocalAssembler::MapToContigs() {
    mapped_f_.resize(contigs_->size());
    mapped_r_.resize(contigs_->size());
    locks_.resize(std::min(contigs_->size(), (size_t)kMaxNumLocks));

    for (auto it = locks_.begin(); it != locks_.end(); ++it) {
        omp_init_lock(&*it);
    }

    max_read_len_ = 1;
    local_range_ = 0;

    for (unsigned lib_id = 0; lib_id < lib_info_.size(); ++lib_id) {
        int local_range = LocalRange_(lib_id);
        bool is_paired = lib_info_[lib_id].is_pe;

        local_range_ = std::max(local_range, local_range_);
        max_read_len_ = std::max(max_read_len_, lib_info_[lib_id].max_read_len);

        MappingRecord rec1, rec2;
        size_t num_added = 0, num_mapped = 0;

        #pragma omp parallel for private(rec1, rec2) reduction(+: num_added, num_mapped) schedule(static, 1)

        for (int64_t i = lib_info_[lib_id].from; i <= lib_info_[lib_id].to; i += 2) {
            bool map1 = MapToHashMapper_(mapper_, i, rec1);
            bool map2 = (i + 1) <= lib_info_[lib_id].to ? MapToHashMapper_(mapper_, i + 1, rec2) : false;

            if (map1) {
                num_added += AddToMappingDeque_(i, rec1, local_range);
                ++num_mapped;

                if (is_paired) {
                    num_added += AddMateToMappingDeque_(i, i + 1, rec1, rec2, map2, local_range);
                }
            }

            if (map2) {
                ++num_mapped;
                num_added += AddToMappingDeque_(i + 1, rec2, local_range);

                if (is_paired) {
                    num_added += AddMateToMappingDeque_(i + 1, i, rec2, rec1, map1, local_range);
                }
            }
        }

        xlog("Lib %d: total %ld reads, aligned %lu, added %lu reads for local assembly\n",
             lib_id, lib_info_[lib_id].to - lib_info_[lib_id].from + 1, num_mapped, num_added);
    }

    for (auto it = locks_.begin(); it != locks_.end(); ++it) {
        omp_destroy_lock(&*it);
    }

    std::vector<omp_lock_t> empty;
    empty.swap(locks_);
}


inline void LaunchIDBA(std::deque<Sequence> &reads, Sequence &contig_end,
                       std::deque<Sequence> &out_contigs,
                       std::deque<ContigInfo> &out_contig_infos,
                       int mink, int maxk, int step) {
    int local_range = contig_end.size();
    HashGraph hash_graph;
    hash_graph.reserve(4 * local_range);

    ContigGraph contig_graph;
    out_contigs.clear();
    out_contig_infos.clear();

    int max_read_len = 0;

    for (unsigned i = 0; i < reads.size(); ++i) {
        max_read_len = std::max(max_read_len, (int)reads[i].size());
    }

    for (int kmer_size = mink; kmer_size <= std::min(maxk, max_read_len); kmer_size += step) {
        int64_t sum = 0;
        hash_graph.set_kmer_size(kmer_size);
        hash_graph.clear();

        for (int64_t i = 0; i < (int64_t)reads.size(); ++i) {
            if ((int)reads[i].size() < kmer_size)
                continue;

            Sequence seq(reads[i]);
            hash_graph.InsertKmers(seq);
            sum += seq.size() - kmer_size + 1;
        }

        Histgram<int> histgram = hash_graph.coverage_histgram();
        double mean = histgram.percentile(1 - 1.0 * local_range / hash_graph.num_vertices());
        double threshold = mean;

        hash_graph.InsertKmers(contig_end);

        for (int64_t i = 0; i < (int64_t)out_contigs.size(); ++i)
            hash_graph.InsertUncountKmers(out_contigs[i]);

        hash_graph.Assemble(out_contigs, out_contig_infos);

        contig_graph.set_kmer_size(kmer_size);
        contig_graph.Initialize(out_contigs, out_contig_infos);
        contig_graph.RemoveDeadEnd(kmer_size * 2);

        contig_graph.RemoveBubble();
        contig_graph.IterateCoverage(kmer_size * 2, 1, threshold);

        contig_graph.Assemble(out_contigs, out_contig_infos);

        if (out_contigs.size() == 1) {
            break;
        }
    }
}

void LocalAssembler::LocalAssemble() {
    output_lock_ = 0;
    int min_num_reads = local_range_ / max_read_len_;

    Sequence seq, contig_end;
    std::deque<Sequence> reads;
    std::deque<Sequence> out_contigs;
    std::deque<ContigInfo> out_contig_infos;

    std::ofstream local_file(local_filename_);
    std::ofstream local_info(local_filename_ + ".info");

    long long num_bases = 0;
    long long num_contigs = 0;

    #pragma omp parallel for private(seq, contig_end, reads, out_contigs, out_contig_infos) schedule(dynamic)

    for (uint64_t cid = 0; cid < contigs_->size(); ++cid) {
        int cl = contigs_->length(cid);

        for (int strand = 0; strand < 2; ++strand) {
            std::deque<uint64_t> &mapped_reads = strand == 0 ? mapped_f_[cid] : mapped_r_[cid];

            if ((int)mapped_reads.size() <= min_num_reads) {
                continue;
            }

            // collect local reads, convert them into Sequence
            reads.clear();

            std::sort(mapped_reads.begin(), mapped_reads.end());
            int last_mapping_pos = -1;
            int pos_count = 0;

            for (uint64_t j = 0; j < mapped_reads.size(); ++j) {
                int pos = mapped_reads[j] >> 49;
                pos_count = pos == last_mapping_pos ? pos_count + 1 : 1;
                last_mapping_pos = pos;

                if (pos_count <= 3) {
                    seq.clear();
                    uint64_t read_id = mapped_reads[j] & ((1ULL << 44) - 1);

                    for (unsigned ri = 0, rsz = reads_->length(read_id); ri < rsz; ++ri) {
                        seq.Append(reads_->get_base(read_id, ri));
                    }

                    reads.push_back(seq);
                }
            }

            contig_end.clear();

            if (strand == 0) {
                for (int j = 0, e = std::min(local_range_, cl); j < e; ++j) {
                    contig_end.Append(contigs_->get_base(cid, j));
                }
            }
            else {
                for (int j = std::max(0,  cl - local_range_); j < cl; ++j) {
                    contig_end.Append(contigs_->get_base(cid, j));
                }
            }

            out_contigs.clear();
            LaunchIDBA(reads, contig_end, out_contigs, out_contig_infos, local_kmin_, local_kmax_, local_step_);

            for (uint64_t j = 0; j < out_contigs.size(); ++j) {
                if (out_contigs[j].size() > (unsigned)min_contig_len_ &&
                        out_contigs[j].size() > (unsigned)local_kmax_) {
                    while (__sync_lock_test_and_set(&output_lock_, 1)) while (output_lock_);

                    WriteFasta(local_file,
                               out_contigs[j],
                               FormatString("lc_%" PRIu64 "_strand_%d_id_%" PRIu64 " flag=0 multi=1",
                                            cid, strand, j));
                    num_contigs++;
                    num_bases += out_contigs[j].size();
                    __sync_lock_release(&output_lock_);
                }
            }
        }
    }

    local_info << num_contigs << ' ' << num_bases << std::endl;

    local_info.close();
    local_file.close();
}