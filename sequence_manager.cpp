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

#include "sequence_manager.h"

#include <assert.h>
#include <zlib.h>
#include <string>

#include "utils.h"
#include "bit_operation.h"

void SequenceManager::set_file(const std::string &file_name) {
    assert(f_type != kMegahitEdges && f_type != kSortedEdges);
    assert(files_.size() == 0);
    assert(kseq_readers_.size() == 0);

    files_.resize(1);
    kseq_readers_.resize(1);

    files_[0] = file_name == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name.c_str(), "r");
    assert(files_[0] != NULL);

    if (f_type == kFastxReads || f_type == kMegahitContigs) {
        kseq_readers_[0] = kseq_init(files_[0]);
        assert(kseq_readers_[0] != NULL);
    }
}

void SequenceManager::set_pe_files(const std::string &file_name1, const std::string &file_name2) {
    assert(f_type == kFastxReads && r_type == kPaired);
    assert(files_.size() == 0);
    assert(kseq_readers_.size() == 0);

    files_.resize(2);
    kseq_readers_.resize(2);

    files_[0] = file_name1 == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name1.c_str(), "r");
    files_[1] = file_name2 == "-" ? gzdopen(fileno(stdin), "r") : gzopen(file_name2.c_str(), "r");

    assert(files_[0] != NULL);

    if (f_type == kFastxReads || f_type == kMegahitContigs) {
        for (int i = 0; i < 2; ++i) {
            kseq_readers_[i] = kseq_init(files_[i]);
            assert(kseq_readers_[i] != NULL);
        }
    }
}

void SequenceManager::set_edge_files(const std::string &file_prefix) {
    assert(f_type == kSortedEdges || f_type == kMegahitEdges);
    assert(files_.size() == 0);
    assert(kseq_readers_.size() == 0);

    assert(!edge_reader_inited_);
    edge_reader_.set_file_prefix(file_prefix);
    edge_reader_.read_info();
    edge_reader_.init_files();

    if (edge_reader_.is_unsorted()) {
        f_type = kMegahitEdges;
    }
    else {
        f_type = kSortedEdges;
    }

    edge_reader_inited_ = true;
}

inline void trimN(const char *s, int len, int &out_bpos, int &out_epos) {
    out_bpos = len;
    out_epos = len;

    int i;

    for (i = 0; i < len; ++i) {
        if (s[i] == 'N' || s[i] == 'n') {
            if (out_bpos < len) {
                break;
            }
        }
        else {
            if (out_bpos == len) {
                out_bpos = i;
            }
        }
    }

    out_epos = i;
}

int64_t SequenceManager::ReadShortReads(int64_t max_num, int64_t max_num_bases, bool append, bool reverse, bool trimN, std::string file_name) {
    if (!append) {
        package_->clear();
    }

    max_num = (max_num + 1) / 2 * 2;
    int64_t num_bases = 0;

    if (f_type == kFastxReads) {
        if (r_type == kPaired) {
            for (int64_t i = 0; i < max_num; i += 2) {
                if (kseq_read(kseq_readers_[0]) >= 0) {
                    if (kseq_read(kseq_readers_[1]) < 0) {
                        xerr_and_exit("File(s) %s: Number of sequences not the same in paired files. Abort.\n", file_name.c_str());
                    }

                    int b0 = 0, e0 = kseq_readers_[0]->seq.l;
                    int b1 = 0, e1 = kseq_readers_[1]->seq.l;

                    if (trimN) {
                        ::trimN(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l, b0, e0);
                        ::trimN(kseq_readers_[1]->seq.s, kseq_readers_[1]->seq.l, b1, e1);
                    }

                    if (reverse) {
                        package_->AppendReverseSeq(kseq_readers_[0]->seq.s + b0, e0 - b0);
                        package_->AppendReverseSeq(kseq_readers_[1]->seq.s + b1, e1 - b1);
                    }
                    else {
                        package_->AppendSeq(kseq_readers_[0]->seq.s + b0, e0 - b0);
                        package_->AppendSeq(kseq_readers_[1]->seq.s + b1, e1 - b1);
                    }

                    num_bases += e0 - b0 + e1 - b1;

                    if (num_bases >= max_num_bases) {
                        return i + 2;
                    }
                }
                else {
                    if (kseq_read(kseq_readers_[1]) >= 0) {
                        xerr_and_exit("File(s) %s: Number of sequences not the same in paired files. Abort.\n", file_name.c_str());
                    }
                    return i;
                }
            }
        }
        else {
            for (int64_t i = 0; i < max_num; ++i) {
                if (kseq_read(kseq_readers_[0]) >= 0) {
                    int b = 0, e = kseq_readers_[0]->seq.l;

                    if (trimN) {
                        ::trimN(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l, b, e);
                    }

                    if (reverse) {
                        package_->AppendReverseSeq(kseq_readers_[0]->seq.s + b, e - b);
                    }
                    else {
                        package_->AppendSeq(kseq_readers_[0]->seq.s + b, e - b);
                    }

                    num_bases += e - b;

                    if (num_bases >= max_num_bases && i % 2 == 1) {
                        return i + 1;
                    }
                }
                else {
                    return i;
                }
            }
        }

        return max_num;
    }
    else if (f_type == kBinaryReads) {
        uint32_t read_len;

        for (int64_t i = 0; i < max_num; ++i) {
            if (gzread(files_[0], &read_len, sizeof(read_len)) == 0) {
                return i;
            }

            int num_words = DivCeiling(read_len, 16);

            if (buf_.size() < (unsigned)num_words) {
                buf_.resize(num_words);
            }

            assert((unsigned)gzread(files_[0], &buf_[0], sizeof(uint32_t) * num_words) == num_words * sizeof(uint32_t));

            if (!reverse) {
                package_->AppendSeq(&buf_[0], read_len);
            }
            else {
                package_->AppendRevSeq(&buf_[0], read_len);
            }

            num_bases += read_len;

            if (read_len >= max_num_bases && i % 2 == 1) {
                return i + 1;
            }
        }

        return max_num;
    }

    assert(false);
}

int64_t SequenceManager::ReadEdges(int64_t max_num, bool append) {
    if (!append) {
        multi_->clear();
        package_->clear();
    }

    if (f_type == kMegahitEdges) {
        for (int64_t i = 0; i < max_num; ++i) {
            uint32_t *next_edge = edge_reader_.NextUnsortedEdge();

            if (next_edge == NULL) {
                return i;
            }

            package_->AppendSeq(next_edge, edge_reader_.kmer_size() + 1);
            multi_->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMulti_t);
        }

        return max_num;
    }
    else if (f_type == kSortedEdges) {
        for (int64_t i = 0; i < max_num; ++i) {
            uint32_t *next_edge = edge_reader_.NextSortedEdge();

            if (next_edge == NULL) {
                return i;
            }

            package_->AppendSeq(next_edge, edge_reader_.kmer_size() + 1);
            multi_->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMulti_t);
        }

        return max_num;
    }

    assert(false);
}


int64_t SequenceManager::ReadEdgesWithFixedLen(int64_t max_num, bool append) {
    if (!append) {
        multi_->clear();
        package_->clear();
    }

    if (f_type == kMegahitEdges) {
        for (int64_t i = 0; i < max_num; ++i) {
            uint32_t *next_edge = edge_reader_.NextUnsortedEdge();

            if (next_edge == NULL) {
                return i;
            }

            package_->AppendFixedLenSeq(next_edge, edge_reader_.kmer_size() + 1);
            multi_->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMulti_t);
        }

        return max_num;
    }
    else if (f_type == kSortedEdges) {
        for (int64_t i = 0; i < max_num; ++i) {
            uint32_t *next_edge = edge_reader_.NextSortedEdge();

            if (next_edge == NULL) {
                return i;
            }

            package_->AppendFixedLenSeq(next_edge, edge_reader_.kmer_size() + 1);
            multi_->push_back(next_edge[edge_reader_.words_per_edge() - 1] & kMaxMulti_t);
        }

        return max_num;
    }

    assert(false);
}

int64_t SequenceManager::ReadMegahitContigs(int64_t max_num, int64_t max_num_bases, bool append, bool reverse,
        int discard_flag, bool extend_loop, bool calc_depth) {
    assert(f_type == kMegahitContigs);
    assert(!(calc_depth && multi_ == NULL));
    assert(!((discard_flag & (contig_flag::kLoop | contig_flag::kIsolated)) && extend_loop)); // loop must be isolated

    if (!append) {
        if (multi_ != NULL) {
            multi_->clear();
        }

        package_->clear();
    }

    int64_t num_bases = 0;

    for (int64_t i = 0; i < max_num; ++i) {
        if (kseq_read(kseq_readers_[0]) >= 0) {
            if ((int)kseq_readers_[0]->seq.l < min_len_ ) {
                --i;
                continue;
            }

            // comment = "flag=x multi=xx.xxxx"
            if (discard_flag & (kseq_readers_[0]->comment.s[5] - '0')) {
                --i;
                continue;
            }

            if (extend_loop && ((kseq_readers_[0]->comment.s[5] - '0') & contig_flag::kLoop)) {
                if (kseq_readers_[0]->seq.l < k_to_ + 1U) {
                    continue;
                }

                std::string ss(kseq_readers_[0]->seq.s);

                for (int i = 0; i < k_to_ - k_from_; ++i) {
                    ss.push_back(ss[i + k_from_]);
                }

                if (reverse) {
                    package_->AppendReverseSeq(ss.c_str(), ss.length());
                }
                else {
                    package_->AppendSeq(ss.c_str(), ss.length());
                }
            }
            else {
                if (reverse) {
                    package_->AppendReverseSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
                }
                else {
                    package_->AppendSeq(kseq_readers_[0]->seq.s, kseq_readers_[0]->seq.l);
                }
            }
            
            float mul = atof(kseq_readers_[0]->comment.s + 13);

            if (multi_ != NULL) {
                multi_->push_back(std::min(kMaxMulti_t, (int)(mul + 0.5)));
            }

            if (float_multi_ != NULL) {
                float_multi_->push_back(mul);
            }

            num_bases += kseq_readers_[0]->seq.l;

            if (num_bases >= max_num_bases) {
                return i + 1;
            }
        }
        else {
            return i;
        }
    }

    return max_num;
}

void SequenceManager::WriteBinarySequences(FILE *file, bool reverse, int64_t from, int64_t to) {
    if (to == -1) {
        to = package_->size() - 1;
    }

    uint32_t len;
    std::vector<uint32_t> s;

    for (int64_t i = from; i <= to; ++i) {
        len = package_->length(i);
        package_->get_seq(s, i);

        if (reverse) {
            for (int j = 0; j < (int)s.size(); ++j) {
                s[j] = bit_operation::Reverse(s[j]);
            }

            for (int j = 0, k = s.size() - 1; j < k; ++j, --k) {
                std::swap(s[j], s[k]);
            }

            int shift = (16 - len % 16)  * 2;

            if (shift != 32) {
                for (int j = 0; j < (int)s.size() - 1; ++j) {
                    s[j] = (s[j] << shift) | (s[j + 1] >> (32 - shift));
                }

                s.back() <<= shift;
            }
        }

        fwrite(&len, sizeof(uint32_t), 1, file);
        fwrite(&s[0], sizeof(uint32_t), s.size(), file);
    }
}

