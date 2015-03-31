#include "cx1_edge2sdbg.h"

#include <omp.h>
#include <string>
#include <vector>

#include "io_utility.h"
#include "utils.h"
#include "packed_reads.h"

#ifndef USE_GPU
#include "lv2_cpu_sort.h"
#else
#include "lv2_gpu_functions.h"
#endif

namespace cx1_edge2sdbg {

// helpers
typedef CX1<edge2sdbg_global_t, kNumBuckets> cx1_t;
typedef CX1<edge2sdbg_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;
typedef CX1<edge2sdbg_global_t, kNumBuckets>::bucketpartition_data_t bucketpartition_data_t;
typedef CX1<edge2sdbg_global_t, kNumBuckets>::outputpartition_data_t outputpartition_data_t;

/**
 * @brief encode edge_id and its offset in one int64_t
 */
inline int64_t EncodeEdgeOffset(int64_t edge_id, int offset, int strand, int k_num_bits) {
    return (edge_id << (k_num_bits + 1)) | (strand << k_num_bits) | offset;
}

/**
 * @brief read reads and search mercy
 */
struct ReadReadsThreadData {
    ReadPackage *read_package;
    gzFile *read_file;
};

static void* ReadReadsThread(void* data) {
    ReadPackage &package = *(((ReadReadsThreadData*)data)->read_package);
    gzFile &read_file = *(((ReadReadsThreadData*)data)->read_file);
    package.ReadBinaryReads(read_file);
    return NULL;
}

/**
 * @brief build lkt for faster binary search for mercy
 */
void InitLookupTable(int64_t *lookup_table, uint32_t *packed_edges, int64_t num_edges, int words_per_edge) {
    memset(lookup_table, 0xFF, sizeof(int64_t) * kLookUpSize * 2);
    uint32_t *edge_p = packed_edges;
    uint32_t cur_prefix = *packed_edges >> kLookUpShift;
    lookup_table[cur_prefix * 2] = 0;
    edge_p += words_per_edge;
    for (int64_t i = 1; i < num_edges; ++i) {
        if ((*edge_p >> kLookUpShift) > cur_prefix) {
            lookup_table[cur_prefix * 2 + 1] = i - 1;
            cur_prefix = (*edge_p >> kLookUpShift);
            lookup_table[cur_prefix * 2] = i;
        } else {
            assert(cur_prefix == (*edge_p >> kLookUpShift));   
        }
        edge_p += words_per_edge;
    }
    lookup_table[cur_prefix * 2 + 1] = num_edges - 1;
}

/**
 * @brief search mercy kmer
 */
inline int64_t BinarySearchKmer(uint32_t *packed_edges, int64_t *lookup_table, int words_per_edge, 
    int words_per_kmer, int last_shift, uint32_t *kmer) {
    // --- first look up ---
    int64_t l = lookup_table[(*kmer >> kLookUpShift) * 2];
    if (l == -1) { return -1; }
    int64_t r = lookup_table[(*kmer >> kLookUpShift) * 2 + 1];
    int64_t mid;

    // --- search the words before the last word ---
    for (int i = 0; i < words_per_kmer - 1; ++i) {
        while (l <= r) {
            mid = (l + r) / 2;
            if (packed_edges[mid * words_per_edge + i] < kmer[i]) {
                l = mid + 1;
            } else if (packed_edges[mid * words_per_edge + i] > kmer[i]) {
                r = mid - 1;
            } else {
                int64_t ll = l, rr = mid, mm;
                while (ll < rr) {
                    mm = (ll + rr) / 2;
                    if (packed_edges[mm * words_per_edge + i] < kmer[i]) {
                        ll = mm + 1;
                    } else {
                        rr = mm;
                    }
                }
                l = ll;

                ll = mid, rr = r;
                while (ll < rr) {
                    mm = (ll + rr + 1) / 2;
                    if (packed_edges[mm * words_per_edge + i] > kmer[i]) {
                        rr = mm - 1;
                    } else {
                        ll = mm;
                    }
                }
                r = rr;
                break;
            }
        }
        if (l > r) { return -1; }
    }

    // --- search the last word ---
    while (l <= r) {
        mid = (l + r) / 2;
        if ((packed_edges[mid * words_per_edge + words_per_kmer - 1] >> last_shift) ==
            (kmer[words_per_kmer - 1] >> last_shift)) {
            return mid;
        } else if ((packed_edges[mid * words_per_edge + words_per_kmer - 1] >> last_shift) > 
            (kmer[words_per_kmer - 1] >> last_shift)) {
            r = mid - 1;
        } else {
            l = mid + 1;
        }
    }
    return -1;
}

/**
 * @brief read candidate reads and search mercy kmer
 * @TODO: many hard-code in this function
 */
void ReadReadsAndGetMercyEdges(edge2sdbg_global_t &globals) {
	int64_t *edge_lookup;
    assert((edge_lookup = (int64_t *) MallocAndCheck(kLookUpSize * 2 * sizeof(int64_t), __FILE__, __LINE__)) != NULL);
    InitLookupTable(edge_lookup, globals.packed_edges, globals.num_edges, globals.words_per_edge);

    uint32_t *packed_edges = globals.packed_edges;
    int64_t *lookup_table = edge_lookup;
    int kmer_k = globals.kmer_k;
    int words_per_edge = globals.words_per_edge;
    const char *edge_file_prefix = globals.input_prefix;

    gzFile candidate_file = gzopen((std::string(globals.input_prefix)+".cand").c_str(), "r");
    ReadPackage read_package[2];
    read_package[0].init(globals.max_read_length);
    read_package[1].init(globals.max_read_length);

    ReadReadsThreadData input_thread_data;
    input_thread_data.read_package = &read_package[0];
    input_thread_data.read_file = &candidate_file;
    int input_thread_idx = 0;
    pthread_t input_thread;

    pthread_create(&input_thread, NULL, ReadReadsThread, &input_thread_data);
    int num_threads = globals.num_cpu_threads - 1;
    omp_set_num_threads(num_threads);

    std::vector<FILE*> out_files;
    for (int i = 0; i < num_threads; ++i) {
        char file_name[10240];
        sprintf(file_name, "%s.mercy.%d", edge_file_prefix, i);
        out_files.push_back(OpenFileAndCheck(file_name, "wb"));
        assert(out_files.back() != NULL);
    }

    // parameters for binary search
    int words_per_kmer = DivCeiling(kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);
    int last_shift_k = (kmer_k * kBitsPerEdgeChar) % kBitsPerEdgeWord;
    if (last_shift_k > 0) {
        last_shift_k = kBitsPerEdgeWord - last_shift_k;
    }
    int words_per_k_plus_one = DivCeiling((kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
    int last_shift_k_plus_one = ((kmer_k + 1) * kBitsPerEdgeChar) % kBitsPerEdgeWord;
    if (last_shift_k_plus_one > 0) {
        last_shift_k_plus_one = kBitsPerEdgeWord - last_shift_k_plus_one;
    }
    // log("%d %d %d %d\n", words_per_kmer, words_per_k_plus_one, last_shift_k, last_shift_k_plus_one);
    uint32_t *kmers = (uint32_t *) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * words_per_edge, __FILE__, __LINE__);
    uint32_t *rev_kmers = (uint32_t *) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * words_per_edge, __FILE__, __LINE__);
    bool *has_ins = (bool*) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * read_package[0].max_read_len, __FILE__, __LINE__);
    bool *has_outs = (bool*) MallocAndCheck(sizeof(uint32_t) * (omp_get_max_threads()) * read_package[0].max_read_len, __FILE__, __LINE__);

    int64_t num_mercy_edges = 0;
    int64_t num_reads = 0;

    while (true) {
        pthread_join(input_thread, NULL);
        ReadPackage &package = read_package[input_thread_idx];
        if (package.num_of_reads == 0) {
            break;
        }

        input_thread_idx ^= 1;
        input_thread_data.read_package = &read_package[input_thread_idx];
        pthread_create(&input_thread, NULL, ReadReadsThread, &input_thread_data);

        num_reads += package.num_of_reads;

#pragma omp parallel for reduction(+:num_mercy_edges)
        for (int read_id = 0; read_id < package.num_of_reads; ++read_id) {
            int read_length = package.length(read_id);
            if (read_length < kmer_k + 2) { continue; }
            bool *has_in = has_ins + omp_get_thread_num() * package.max_read_len;
            bool *has_out = has_outs + omp_get_thread_num() * package.max_read_len;
            memset(has_in, 0, sizeof(bool) * (read_length - kmer_k + 1));
            memset(has_out, 0, sizeof(bool) * (read_length - kmer_k + 1));
            // construct the first kmer
            uint32_t *kmer = kmers + words_per_edge * omp_get_thread_num();
            uint32_t *rev_kmer = rev_kmers + words_per_edge * omp_get_thread_num();
            memcpy(kmer, package.GetReadPtr(read_id), sizeof(uint32_t) * words_per_k_plus_one);
            // construct the rev_kmer
            for (int i = 0; i < words_per_kmer; ++i) {
                rev_kmer[words_per_kmer - 1 - i] = ~ mirror(kmer[i]);
            }
            for (int i = 0; i < words_per_kmer; ++i) {
                rev_kmer[i] <<= last_shift_k;
                rev_kmer[i] |= (i == words_per_kmer - 1) ? 0 : (rev_kmer[i + 1] >> (kBitsPerEdgeWord - last_shift_k));
            }

            int last_index = std::min(read_length - 1, kCharsPerEdgeWord * words_per_k_plus_one - 1);

            // first determine which kmer has in or out
            for (int first_index = 0; first_index + kmer_k <= read_length; ++first_index) {
                if (!has_in[first_index]) {
                    // search the reverse complement
                    if (BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                            words_per_kmer, last_shift_k, rev_kmer) != -1) {
                        has_in[first_index] = true;
                    } else {
                        // check whether it has incomings
                        int last_char = kmer[words_per_k_plus_one - 1] & 3;
                        for (int i = words_per_k_plus_one - 1; i > 0; --i) {
                            kmer[i] = (kmer[i] >> 2) | (kmer[i - 1] << 30);
                        }
                        kmer[0] >>= 2;
                        // set the highest char to c
                        for (int c = 0; c < 4; ++c) {
                            kmer[0] &= 0x3FFFFFFF;
                            kmer[0] |= c << 30;
                            if (kmer[0] > rev_kmer[0]) {
                                break;
                            }
                            if (BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                     words_per_k_plus_one, last_shift_k_plus_one, kmer) != -1) {
                                has_in[first_index] = true;
                                break;
                            }
                        }
                        for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                            kmer[i] = (kmer[i] << 2) | (kmer[i + 1] >> 30);
                        }
                        kmer[words_per_k_plus_one - 1] = (kmer[words_per_k_plus_one - 1] << 2) | last_char;
                    }
                }

                if (true) {
                    // check whether it has outgoing
                    int64_t search_idx = BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                                          words_per_kmer, last_shift_k, kmer);
                    if (search_idx != -1) {
                        has_out[first_index] = true;
                        // a quick check whether next has in
                        if (first_index + kmer_k < read_length && 
                            (packed_edges[search_idx * words_per_edge + words_per_k_plus_one - 1] >> last_shift_k_plus_one) ==
                            (kmer[words_per_k_plus_one - 1] >> last_shift_k_plus_one)) {
                            has_in[first_index + 1] = true;
                        }
                    } else {
                        // search the rc
                        int rc_last_char = rev_kmer[words_per_k_plus_one - 1] & 3;
                        for (int i = words_per_k_plus_one - 1; i > 0; --i) {
                            rev_kmer[i] = (rev_kmer[i] >> 2) | (rev_kmer[i - 1] << 30);
                        }
                        rev_kmer[0] >>= 2;
                        int next_c = first_index + kmer_k < read_length ?
                                     (3 - package.CharAt(read_id, first_index + kmer_k)) :
                                     3;
                        rev_kmer[0] &= 0x3FFFFFFF;
                        rev_kmer[0] |= next_c << 30;
                        if (rev_kmer[0] <= kmer[0] &&
                            BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                words_per_k_plus_one, last_shift_k_plus_one, rev_kmer) != -1) {
                            has_out[first_index] = true;
                            has_in[first_index + 1] = true;
                        }

                        for (int c = 0; !has_out[first_index] && c < 4; ++c) {
                            if (c == next_c) { continue; }
                            rev_kmer[0] &= 0x3FFFFFFF;
                            rev_kmer[0] |= c << 30;
                            if (rev_kmer[0] > kmer[0]) {
                                break;
                            }
                            if (BinarySearchKmer(packed_edges, lookup_table, words_per_edge, 
                                    words_per_k_plus_one, last_shift_k_plus_one, rev_kmer) != -1) {
                                has_out[first_index] = true;
                                break;
                            }
                        }
                        for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                            rev_kmer[i] = (rev_kmer[i] << 2) | (rev_kmer[i + 1] >> 30);
                        }
                        rev_kmer[words_per_k_plus_one - 1] = (rev_kmer[words_per_k_plus_one - 1] << 2) | rc_last_char;
                    }
                }

                // shift kmer and rev_kmer
                for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                    kmer[i] = (kmer[i] << 2) | (kmer[i + 1] >> 30);
                }
                kmer[words_per_k_plus_one - 1] <<= 2;
                if (++last_index < read_length) {
                    kmer[words_per_k_plus_one - 1] |= package.CharAt(read_id, last_index);
                }

                for (int i = words_per_k_plus_one - 1; i > 0; --i) {
                    rev_kmer[i] = (rev_kmer[i] >> 2) | (rev_kmer[i - 1] << 30);
                }
                rev_kmer[0] = (rev_kmer[0] >> 2) | ((3 - package.CharAt(read_id, first_index + kmer_k)) << 30);
            }

            // adding mercy edges
            int last_no_out = -1;
            std::vector<bool> is_mercy_edges(read_length - kmer_k, false);
            for (int i = 0; i + kmer_k <= read_length; ++i) {
                switch (has_in[i] | (int(has_out[i]) << 1)) {
                    case 1: { // has incoming only
                        last_no_out = i;
                        break;
                    }
                    case 2: { // has outgoing only
                        if (last_no_out >= 0) {
                            for (int j = last_no_out; j < i; ++j) {
                                is_mercy_edges[j] = true;
                            }
                            num_mercy_edges += i - last_no_out;
                        }
                        last_no_out = -1;
                        break;
                    }
                    case 3: { // has in and out
                        last_no_out = -1;
                        break;
                    }
                    default: {
                        // do nothing
                        break;
                    }
                }
            }
            
            memcpy(kmer, package.GetReadPtr(read_id), sizeof(uint32_t) * words_per_k_plus_one);
            last_index = std::min(read_length - 1, kCharsPerEdgeWord * words_per_k_plus_one - 1);
            for (int i = 0; i + kmer_k < read_length; ++i) {
                if (is_mercy_edges[i]) {
                    uint32_t last_word = kmer[words_per_k_plus_one - 1];
                    kmer[words_per_k_plus_one - 1] >>= last_shift_k_plus_one;
                    kmer[words_per_k_plus_one - 1] <<= last_shift_k_plus_one;
                    for (int j = words_per_k_plus_one; j < words_per_edge; ++j) {
                        kmer[j] = 0;
                    }
                    if (globals.mult_mem_type == 0) {
                        kmer[words_per_edge - 1] |= 1; // WARNING: only accurate when m=2, but I think doesn't matter a lot
                    }
                    fwrite(kmer, sizeof(uint32_t), words_per_edge, out_files[omp_get_thread_num()]);
                    if (globals.mult_mem_type > 0) {
                        uint32_t kMercyMult = 1;
                        fwrite(&kMercyMult, sizeof(uint32_t), 1, out_files[omp_get_thread_num()]);
                    }
                    kmer[words_per_k_plus_one - 1] = last_word;
                }

                for (int i = 0; i < words_per_k_plus_one - 1; ++i) {
                    kmer[i] = (kmer[i] << 2) | (kmer[i + 1] >> 30);
                }
                kmer[words_per_k_plus_one - 1] <<= 2;
                if (++last_index < read_length) {
                    kmer[words_per_k_plus_one - 1] |= package.CharAt(read_id, last_index);
                }
            }
        }
        if (num_reads % (16 * package.kMaxNumReads) == 0) {
            if (cx1_t::kCX1Verbose >= 4) {
                log("[B::%s] Number of reads: %ld, Number of mercy edges: %ld\n", __func__, num_reads, num_mercy_edges);
            }
        }
    }

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] Number of reads: %ld, Number of mercy edges: %ld\n", __func__, num_reads, num_mercy_edges);
    }

    free(kmers);
    free(rev_kmers);
    free(has_ins);
    free(has_outs);
    free(edge_lookup);
    for (unsigned i = 0; i < out_files.size(); ++i) {
        fclose(out_files[i]);
    }
}

/**
 * @brief read mercy edges from disk
 */
int64_t ReadMercyEdges(edge2sdbg_global_t &globals) {
    // --- init reader ---
    EdgeReader edge_reader;
    edge_reader.InitUnsorted((std::string(globals.input_prefix) + ".mercy").c_str(),
                             globals.num_cpu_threads - 1,
                             globals.kmer_k,
                             globals.words_per_edge + (globals.mult_mem_type > 0));

    // --- read mercy edges ---
    edge_word_t *edge_p = globals.packed_edges + globals.num_edges * globals.words_per_edge;
    int64_t bytes_per_edge = globals.words_per_edge * sizeof(edge_word_t) + globals.mult_mem_type;
    int64_t max_num_edges = globals.host_mem * 0.9 / bytes_per_edge; // TODO: more accurate
    int64_t num_edges = globals.num_edges;
    while (edge_reader.NextEdgeUnsorted(edge_p)) {
        if (num_edges >= globals.capacity) {
            if (num_edges >= max_num_edges) {
                err("[B::%s ERROR] reach max_num_edges: %ld... No enough memory to build the graph.\n", __func__, max_num_edges);
                num_edges = globals.num_edges; // reset
                exit(1);
            }

            globals.capacity = std::min(max_num_edges, globals.capacity * 2);
            edge_word_t *new_packed_edge = (edge_word_t*) realloc(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity);
            if (new_packed_edge == NULL) {
                err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                    __func__, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity, globals.host_mem);
                exit(1);
            }
            globals.packed_edges = new_packed_edge;
            edge_p = new_packed_edge + num_edges * globals.words_per_edge;

            if (globals.mult_mem_type == 1) {
                uint8_t *new_multi8 = (uint8_t*) realloc(globals.multiplicity8, sizeof(uint8_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint8_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity8 = new_multi8;
            } else if (globals.mult_mem_type == 2) {
                uint16_t *new_multi16 = (uint16_t*) realloc(globals.multiplicity16, sizeof(uint16_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint16_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity16 = new_multi16;
            }
        }

        ++num_edges;
        if (globals.mult_mem_type == 1) {
            edge_p[globals.words_per_edge - 1] &= 0xFFFFFF00U;
            edge_p[globals.words_per_edge - 1] |= (edge_p[globals.words_per_edge] >> 8) & 0xFFU;
            globals.multiplicity8[num_edges - 1] = edge_p[globals.words_per_edge] & 0xFFU;
        } else if (globals.mult_mem_type == 2) {
            globals.multiplicity16[num_edges - 1] = edge_p[globals.words_per_edge] & kMaxMulti_t;
        }

        edge_p += globals.words_per_edge;
    }
    edge_reader.destroy();

    // --- realloc ---
    globals.packed_edges = (edge_word_t *) ReAllocAndCheck(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * num_edges, __FILE__, __LINE__);
    if (globals.mult_mem_type == 1) {
        globals.multiplicity8 = (uint8_t*) ReAllocAndCheck(globals.multiplicity8, sizeof(uint8_t) * num_edges, __FILE__, __LINE__);
    } else if (globals.mult_mem_type == 2) {
        globals.multiplicity16 = (uint16_t*) ReAllocAndCheck(globals.multiplicity16, sizeof(uint16_t) * num_edges, __FILE__, __LINE__);
    }


    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] Number of mercy edges: %lld\n", __func__, num_edges - globals.num_edges);
    }
    globals.num_edges = num_edges;
    globals.mem_packed_edges = bytes_per_edge * num_edges;
    return num_edges;
}

inline bool IsDiffKMinusOneMer(edge_word_t *item1, edge_word_t *item2, int64_t spacing, int kmer_k) {
    // mask extra bits
    int chars_in_last_word = (kmer_k - 1) % kCharsPerEdgeWord;
    int num_full_words = (kmer_k - 1) / kCharsPerEdgeWord;
    if (chars_in_last_word > 0) {
        edge_word_t w1 = item1[num_full_words * spacing];
        edge_word_t w2 = item2[num_full_words * spacing];
        if ((w1 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar) != (w2 >> (kCharsPerEdgeWord - chars_in_last_word) * kBitsPerEdgeChar)) {
            return true;
        }
    } 

    for (int i = num_full_words - 1; i >= 0; --i) {
        if (item1[i * spacing] != item2[i * spacing]) {
            return true;
        }
    }
    return false;
}

inline int ExtractFirstChar(edge_word_t *item) {
    return *item >> kTopCharShift;
}

inline int Extract_a(edge_word_t *item, int num_words, int64_t spacing, int kmer_k) {
    int non_dollar = (item[(num_words - 1) * spacing] >> (kBWTCharNumBits + kBitsPerMulti_t)) & 1;
    if (non_dollar) {
        int which_word = (kmer_k - 1) / kCharsPerEdgeWord;
        int word_index = (kmer_k - 1) % kCharsPerEdgeWord;
        return (item[which_word * spacing] >> (kCharsPerEdgeWord - 1 - word_index) * kBitsPerEdgeChar) & kEdgeCharMask;
    } else {
        return kSentinelValue;
    }
}

inline int Extract_b(edge_word_t *item, int num_words, int64_t spacing) {
    return (item[(num_words - 1) * spacing] >> kBitsPerMulti_t) & ((1 << kBWTCharNumBits) - 1);
}

inline int ExtractCounting(edge_word_t *item, int num_words, int64_t spacing) {
    return item[(num_words - 1) * spacing] & kMaxMulti_t; 
}

// cx1 core functions
int64_t encode_lv1_diff_base(int64_t read_id, edge2sdbg_global_t &g) {
	return EncodeEdgeOffset(read_id, 0, 0, g.k_num_bits);
}

void read_edges_and_prepare(edge2sdbg_global_t &globals) {
    // --- init reader ---
    EdgeReader edge_reader;
    edge_reader.init((std::string(globals.input_prefix) + ".edges").c_str(), globals.num_edge_files);
    // --- calc memory type and max_num_edge ---
    globals.kmer_k = edge_reader.kmer_k;
    globals.words_per_edge = DivCeiling(globals.kmer_k + 1, kCharsPerEdgeWord);
    int free_bits_in_edge = globals.words_per_edge * kBitsPerEdgeWord - (globals.kmer_k + 1) * kBitsPerEdgeChar;
    if (free_bits_in_edge >= kBitsPerMulti_t) {
        globals.mult_mem_type = 0;
    } else if (free_bits_in_edge >= 8) {
        globals.mult_mem_type = 1;
    } else {
        globals.mult_mem_type = 2;
    }
    int64_t bytes_per_edge = globals.words_per_edge * sizeof(edge_word_t) + globals.mult_mem_type;
    int64_t max_num_edges = globals.host_mem * 0.9 / bytes_per_edge; // TODO: more accurate

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] kmer_k: %d, words_per_edge: %d\n", __func__, globals.kmer_k, globals.words_per_edge);
        log("[B::%s] Max host mem: %ld, max number of edges can be loaded: %lld\n", __func__, globals.host_mem, (long long)max_num_edges);
    }

    // --- alloc memory for edges ---
    globals.capacity = std::min(max_num_edges, int64_t(10485760)); // 10M
    globals.packed_edges = (edge_word_t *) MallocAndCheck(sizeof(edge_word_t) * globals.words_per_edge * globals.capacity, __FILE__, __LINE__);
    if (globals.mult_mem_type == 1) {
        globals.multiplicity8 = (uint8_t*) MallocAndCheck(sizeof(uint8_t) * globals.capacity, __FILE__, __LINE__);
    } else if (globals.mult_mem_type == 2) {
        globals.multiplicity16 = (uint16_t*) MallocAndCheck(sizeof(uint16_t) * globals.capacity, __FILE__, __LINE__);
    }

    // --- read edges ---
    edge_word_t *edge_p = globals.packed_edges;
    int64_t num_edges = 0;
    while (edge_reader.NextEdge(edge_p)) {
        if (num_edges >= globals.capacity) {
            if (num_edges >= max_num_edges) {
                err("[B::%s ERROR] reach max_num_edges: %ld... No enough memory to build the graph.\n", __func__, max_num_edges);
                exit(1);
            }

            globals.capacity = std::min(max_num_edges, globals.capacity * 2);
            edge_word_t *new_packed_edge = (edge_word_t*) realloc(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity);
            if (new_packed_edge == NULL) {
                err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                    __func__, sizeof(edge_word_t) * globals.words_per_edge * globals.capacity, globals.host_mem);
                exit(1);
            }
            globals.packed_edges = new_packed_edge;
            edge_p = new_packed_edge + num_edges * globals.words_per_edge;

            if (globals.mult_mem_type == 1) {
                uint8_t *new_multi8 = (uint8_t*) realloc(globals.multiplicity8, sizeof(uint8_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint8_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity8 = new_multi8;
            } else if (globals.mult_mem_type == 2) {
                uint16_t *new_multi16 = (uint16_t*) realloc(globals.multiplicity16, sizeof(uint16_t) * globals.capacity);
                if (new_packed_edge == NULL) {
                    err("[B::%s ERROR] Fail to alloc %lld bytes; the setting of host_mem = %lld is too big?\n", 
                        __func__, sizeof(edge_word_t) * sizeof(uint16_t) * globals.capacity, globals.host_mem);
                    exit(1);
                }
                globals.multiplicity16 = new_multi16;
            }
        }

        ++num_edges;
        if (globals.mult_mem_type == 1) {
            edge_p[globals.words_per_edge - 1] &= 0xFFFFFF00U;
            edge_p[globals.words_per_edge - 1] |= (edge_p[globals.words_per_edge] >> 8) & 0xFFU;
            globals.multiplicity8[num_edges - 1] = edge_p[globals.words_per_edge] & 0xFFU;
        } else if (globals.mult_mem_type == 2) {
            globals.multiplicity16[num_edges - 1] = edge_p[globals.words_per_edge] & kMaxMulti_t;
        }

        edge_p += globals.words_per_edge;
    }
    edge_reader.destroy();

    // --- realloc if no mercy ---
    if (!globals.need_mercy) {
        globals.packed_edges = (edge_word_t *) ReAllocAndCheck(globals.packed_edges, sizeof(edge_word_t) * globals.words_per_edge * num_edges, __FILE__, __LINE__);
        if (globals.mult_mem_type == 1) {
            globals.multiplicity8 = (uint8_t*) ReAllocAndCheck(globals.multiplicity8, sizeof(uint8_t) * num_edges, __FILE__, __LINE__);
        } else if (globals.mult_mem_type == 2) {
            globals.multiplicity16 = (uint16_t*) ReAllocAndCheck(globals.multiplicity16, sizeof(uint16_t) * num_edges, __FILE__, __LINE__);
        }
        globals.mem_packed_edges = bytes_per_edge * num_edges;
    }
    globals.num_edges = num_edges;

    if (cx1_t::kCX1Verbose >= 3) {
        log("[B::%s] Number of solid edges: %lld\n", __func__, num_edges);
    }

    // --- mercy ---
    if (globals.need_mercy) {
    	xtimer_t timer;
        if (cx1_t::kCX1Verbose >= 3) {
	    	timer.reset();
	        timer.start();
            log("[B::%s] Adding mercy edges...\n", __func__);
        }
        ReadReadsAndGetMercyEdges(globals);
        ReadMercyEdges(globals);

        if (cx1_t::kCX1Verbose >= 3) {
        	timer.stop();
            log("[B::%s] Done. Time elapsed: %.4lf\n", __func__, timer.elapsed());
        }
    }

    // --- compute k_num_bits ---
    {
        globals.k_num_bits = 0;
        int len = 1;
        while (len < globals.kmer_k + 1) {
            globals.k_num_bits++;
            len *= 2;
        }
    }

    // --- set cx1 param ---
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_output_threads_ = globals.num_output_threads;
    globals.cx1.num_items_ = globals.num_edges;
}

void* lv0_calc_bucket_size(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
    edge2sdbg_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, kNumBuckets * sizeof(int64_t));
    edge_word_t *edge_p = globals.packed_edges + rp.rp_start_id * globals.words_per_edge;
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, edge_p += globals.words_per_edge) {
        edge_word_t key = 0; // $$$$$$$$
        edge_word_t *word_p = edge_p;
        edge_word_t w = *(word_p++);
        // build initial partial key
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + (w >> kTopCharShift) + 1;
            w <<= kBitsPerEdgeChar;
        }
        // 3 edges Sb$, aSb, $aS
        for (int i = kBucketPrefixLength - 1; i <= kBucketPrefixLength + 1; ++i) {
            if (i % kCharsPerEdgeWord == 0) {
                w = *(word_p++);
            }
            key = (key * kBucketBase + (w >> kTopCharShift) + 1) % kNumBuckets;
            w <<= kBitsPerEdgeChar;
            bucket_sizes[key]++;
        }

        // reverse complement very sucking
        key = 0;
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1; // complement
        }
        for (int i = kBucketPrefixLength - 1; i <= kBucketPrefixLength + 1; ++i) {
            key = key * kBucketBase + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1;
            key %= kNumBuckets;
            bucket_sizes[key]++;
        }
    }
    return NULL;
}

void init_global_and_set_cx1(edge2sdbg_global_t &globals) {
    // --- calculate lv2 memory ---
    globals.max_bucket_size = *std::max_element(globals.cx1.bucket_sizes_, globals.cx1.bucket_sizes_ + kNumBuckets);
    globals.tot_bucket_size = 0;
    for (int i = 0; i < kNumBuckets; ++i) { globals.tot_bucket_size += globals.cx1.bucket_sizes_[i]; }
#ifndef USE_GPU
    globals.cx1.max_lv2_items_ = std::max(globals.max_bucket_size, kMinLv2BatchSize);
#else
    int64_t lv2_mem = globals.gpu_mem - 1073741824; // should reserver ~1G for GPU sorting
    globals.cx1.max_lv2_items_ = std::min(lv2_mem / cx1_t::kGPUBytePerItem, std::max(globals.max_bucket_size, kMinLv2BatchSizeGPU));
    if (globals.max_bucket_size > globals.cx1.max_lv2_items_) {
        err("[ERROR B::%s] Bucket too large for GPU: contains %lld items. Please try CPU version.\n", __func__, globals.max_bucket_size);
        // TODO: auto switch to CPU version
        exit(1);
    }
#endif
    globals.words_per_substring = DivCeiling(globals.kmer_k * kBitsPerEdgeChar + kBWTCharNumBits + 1 + kBitsPerMulti_t, kBitsPerEdgeWord);
    globals.words_per_dummy_node = DivCeiling(globals.kmer_k * kBitsPerEdgeChar, kBitsPerEdgeWord);
    // lv2 bytes: substring (double buffer), permutation, aux
    int64_t lv2_bytes_per_item = (globals.words_per_substring * sizeof(edge_word_t) + sizeof(uint32_t)) * 2 + sizeof(unsigned char);
#ifndef USE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // simulate GPU
#endif

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] %d words per substring, k_num_bits: %d, words per dummy node ($v): %d\n", __func__, globals.words_per_substring, globals.k_num_bits, globals.words_per_dummy_node);
    }

    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem
                         - globals.mem_packed_edges
                         - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1);
    int64_t min_lv1_items = globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5);
    int64_t min_lv2_items = std::max(globals.max_bucket_size, kMinLv2BatchSize);

    if (globals.mem_flag == 1) {
        // auto set memory
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv2_items_, int64_t(globals.tot_bucket_size / (kDefaultLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.cx1.max_lv2_items_ * lv2_bytes_per_item;
        if (mem_needed > mem_remained) {
        	globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        }
    } else if (globals.mem_flag == 0) {
        // min memory
        globals.cx1.max_lv1_items_ = std::max(globals.cx1.max_lv2_items_, int64_t(globals.tot_bucket_size / (kMaxLv1ScanTime - 0.5)));
        int64_t mem_needed = globals.cx1.max_lv1_items_ * cx1_t::kLv1BytePerItem + globals.cx1.max_lv2_items_ * lv2_bytes_per_item;
        if (mem_needed > mem_remained) {
        	globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        } else {
        	globals.cx1.adjust_mem(mem_needed, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
        }
    } else {
        // use all
        globals.cx1.adjust_mem(mem_remained, lv2_bytes_per_item, min_lv1_items, min_lv2_items);
    }

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] Memory for edges: %lld\n", __func__, globals.mem_packed_edges);
        log("[B::%s] max # lv.1 items = %lld\n", __func__, globals.cx1.max_lv1_items_);
        log("[B::%s] max # lv.2 items = %lld\n", __func__, globals.cx1.max_lv2_items_);
    }

    // --- alloc memory ---
    globals.lv1_items = (int*) MallocAndCheck(globals.cx1.max_lv1_items_ * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (edge_word_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_substrings_db = (edge_word_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_aux = (unsigned char*) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(unsigned char), __FILE__, __LINE__);
 #ifndef USE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.cx1.max_lv2_items_, __FILE__, __LINE__);
 #else
    alloc_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2, (size_t)globals.cx1.max_lv2_items_);
 #endif

	// --- init lock ---
    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL);

    // --- init stat ---
    globals.cur_prefix = -1;
    globals.cur_suffix_first_char = -1;
    globals.num_ones_in_last = 0;
    globals.total_number_edges = 0;
    globals.num_dollar_nodes = 0;
    memset(globals.num_chars_in_w, 0, sizeof(globals.num_chars_in_w));

	// --- init output ---
    globals.sdbg_writer.init((std::string(globals.output_prefix)+".w").c_str(),
        (std::string(globals.output_prefix)+".last").c_str(),
        (std::string(globals.output_prefix)+".isd").c_str());
    globals.dummy_nodes_writer.init((std::string(globals.output_prefix)+".dn").c_str());
    globals.output_f_file = OpenFileAndCheck((std::string(globals.output_prefix)+".f").c_str(), "w");
    globals.output_multiplicity_file = OpenFileAndCheck((std::string(globals.output_prefix)+".mul").c_str(), "wb");
    globals.output_multiplicity_file2 = OpenFileAndCheck((std::string(globals.output_prefix)+".mul2").c_str(), "wb");
    // --- write header ---
    fprintf(globals.output_f_file, "-1\n");
    globals.dummy_nodes_writer.output(globals.words_per_dummy_node);
}

void* lv1_fill_offset(void* _data) {
    readpartition_data_t &rp = *((readpartition_data_t*) _data);
    edge2sdbg_global_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *)MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);
    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;
    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    edge_word_t *edge_p = globals.packed_edges + rp.rp_start_id * globals.words_per_edge;

    // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)                                   \
    do {                                                                \
      if (((key - globals.cx1.lv1_start_bucket_) ^ (key - globals.cx1.lv1_end_bucket_)) & kSignBitMask) { \
        int64_t full_offset = EncodeEdgeOffset(edge_id, offset, strand, globals.k_num_bits); \
        int64_t differential = full_offset - prev_full_offsets[key];      \
        if (differential > cx1_t::kDifferentialLimit) {                      \
          pthread_mutex_lock(&globals.lv1_items_scanning_lock); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = -globals.cx1.lv1_items_special_.size() - 1; \
          globals.cx1.lv1_items_special_.push_back(full_offset);                  \
          pthread_mutex_unlock(&globals.lv1_items_scanning_lock); \
        } else {                                                              \
          assert(differential >= 0); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = (int) differential; \
        } \
        prev_full_offsets[key] = full_offset;                           \
      }                                                                 \
    } while (0)
    // ^^^^^ why is the macro surrounded by a do-while? please ask Google
    // =========== end macro ==========================

    for (int64_t edge_id = rp.rp_start_id; edge_id < rp.rp_end_id; ++edge_id, edge_p += globals.words_per_edge) {
        edge_word_t key = 0; // $$$$$$$$
        edge_word_t *word_p = edge_p;
        edge_word_t w = *(word_p++);
        // build initial partial key
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + (w >> kTopCharShift) + 1;
            w <<= kBitsPerEdgeChar;
        }
        // 3 edges Sb$, aSb, $aS
        for (int i = kBucketPrefixLength - 1; i <= kBucketPrefixLength + 1; ++i) {
            if (i % kCharsPerEdgeWord == 0) {
                w = *(word_p++);
            }
            key = (key * kBucketBase + (w >> kTopCharShift) + 1) % kNumBuckets;
            w <<= kBitsPerEdgeChar;
            CHECK_AND_SAVE_OFFSET(i - kBucketPrefixLength + 1, 0);
        }

        // reverse complement very sucking
        key = 0;
        for (int i = 0; i < kBucketPrefixLength - 1; ++i) {
            key = key * kBucketBase + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1; // complement
        }
        for (int i = kBucketPrefixLength - 1; i <= kBucketPrefixLength + 1; ++i) {
            key = (key * kBucketBase) + (3 - ExtractNthChar(edge_p, globals.kmer_k - i)) + 1;
            key %= kNumBuckets;
            CHECK_AND_SAVE_OFFSET(i - kBucketPrefixLength + 1, 1);
        }
    }
#undef CHECK_AND_SAVE_OFFSET

    free(prev_full_offsets);
    return NULL;
}

void* lv2_extract_substr(void* _data) {
    bucketpartition_data_t &bp = *((bucketpartition_data_t*) _data);
    edge2sdbg_global_t &globals = *(bp.globals);
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[ bp.bp_start_bucket ];
    int64_t offset_mask = (1 << globals.k_num_bits) - 1; // 0000....00011..11
    edge_word_t *substrings_p = globals.lv2_substrings +
                         (globals.cx1.rp_[0].rp_bucket_offsets[ bp.bp_start_bucket ] - globals.cx1.rp_[0].rp_bucket_offsets[ globals.cx1.lv2_start_bucket_ ]);
    for (int bucket = bp.bp_start_bucket; bucket < bp.bp_end_bucket; ++bucket) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
            int num = globals.cx1.rp_[t].rp_bucket_sizes[bucket];
            for (int i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                } else {
                    full_offset = globals.cx1.lv1_items_special_[-1 - *(lv1_p++)];
                }
                int64_t edge_id = full_offset >> (1 + globals.k_num_bits);
                int offset = full_offset & offset_mask;
                int strand = (full_offset >> globals.k_num_bits) & 1;
                edge_word_t *edge_p = globals.packed_edges + edge_id * globals.words_per_edge;
                int num_chars_to_copy = globals.kmer_k - (offset >= 2);
                int counting = 0;
                if (offset == 1) {
                    switch (globals.mult_mem_type) {
                      case 0:
                        counting = *(edge_p + globals.words_per_edge - 1) & kMaxMulti_t;
                        break;
                      case 1:
                        counting = ((*(edge_p + globals.words_per_edge - 1) & 0xFF) << 8) | globals.multiplicity8[edge_id];
                        break;
                      case 2:
                        counting = globals.multiplicity16[edge_id];
                        break;
                      default: assert(false);
                    }
                }
                if (strand == 0) {
                    // copy counting and W char
                    int prev_char;
                    if (offset == 0) {
                        assert(num_chars_to_copy == globals.kmer_k);
                        prev_char = kSentinelValue;
                    } else {
                        prev_char = ExtractNthChar(edge_p, offset - 1);
                    }
                    
                    CopySubstring(substrings_p, edge_p, offset, num_chars_to_copy,
                                  globals.cx1.lv2_num_items_, globals.words_per_edge, globals.words_per_substring);

                    edge_word_t *last_word = substrings_p + int64_t(globals.words_per_substring - 1) * globals.cx1.lv2_num_items_;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << (kBWTCharNumBits + kBitsPerMulti_t);
                    *last_word |= prev_char << kBitsPerMulti_t;
                    *last_word |= std::min(counting, kMaxMulti_t);
                } else {
                    int prev_char;
                    if (offset == 0) {
                        assert(num_chars_to_copy == globals.kmer_k);
                        prev_char = kSentinelValue;
                    } else {
                        prev_char = 3 - ExtractNthChar(edge_p, globals.kmer_k - (offset - 1));
                    }

                    offset = offset == 0 ? 1 : 0; // convert to normal
                    CopySubstringRC(substrings_p, edge_p, offset, num_chars_to_copy, 
                                    globals.cx1.lv2_num_items_, globals.words_per_edge, globals.words_per_substring);

                    edge_word_t *last_word = substrings_p + int64_t(globals.words_per_substring - 1) * globals.cx1.lv2_num_items_;
                    *last_word |= int(num_chars_to_copy == globals.kmer_k) << (kBWTCharNumBits + kBitsPerMulti_t);
                    *last_word |= prev_char << kBitsPerMulti_t;
                    *last_word |= std::min(counting, kMaxMulti_t);
                }
                substrings_p++;
            }
        }
    }
    return NULL;
}

void lv2_sort(edge2sdbg_global_t &globals) {
	xtimer_t local_timer;
#ifndef USE_GPU
    if (cx1_t::kCX1Verbose >= 4) {
	    local_timer.reset();
	    local_timer.start();
    }
    omp_set_num_threads(globals.num_cpu_threads - globals.num_output_threads);
    lv2_cpu_sort(globals.lv2_substrings, globals.permutation, globals.cpu_sort_space, globals.words_per_substring, globals.cx1.lv2_num_items_);
    omp_set_num_threads(globals.num_cpu_threads);
    local_timer.stop();

    if (cx1_t::kCX1Verbose >= 4) {
        log("[B::%s] Sorting substrings with CPU...done. Time elapsed: %.4lf\n", __func__, local_timer.elapsed());
    }
#else
    if (cx1_t::kCX1Verbose >= 4) {
	    local_timer.reset();
	    local_timer.start();
    }
    lv2_gpu_sort(globals.lv2_substrings, globals.permutation, globals.words_per_substring, globals.cx1.lv2_num_items_,
                 globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);

    if (cx1_t::kCX1Verbose >= 4) {
    	local_timer.stop();
        log("[B::%s] Sorting substrings with GPU...done. Time elapsed: %.4lf\n", __func__, local_timer.elapsed());
    }
#endif
}

void lv2_pre_output_partition(edge2sdbg_global_t &globals) {
	// swap double buffers
    globals.lv2_num_items_db = globals.cx1.lv2_num_items_;
    std::swap(globals.lv2_substrings_db, globals.lv2_substrings);
    std::swap(globals.permutation_db, globals.permutation);

	// distribute threads
    int64_t last_end_index = 0;
    int64_t items_per_thread = globals.lv2_num_items_db / globals.num_output_threads;

    for (int t = 0; t < globals.num_output_threads - 1; ++t) {
        int64_t this_start_index = last_end_index;
        int64_t this_end_index = this_start_index + items_per_thread;
        if (this_end_index > globals.lv2_num_items_db) { this_end_index = globals.lv2_num_items_db; }
        if (this_end_index > 0) {
            while (this_end_index < globals.lv2_num_items_db) {
                edge_word_t *prev_item = globals.lv2_substrings_db + globals.permutation_db[this_end_index - 1];
                edge_word_t *item = globals.lv2_substrings_db + globals.permutation_db[this_end_index];
                if (IsDiffKMinusOneMer(item, prev_item, globals.lv2_num_items_db, globals.kmer_k)) {
                    break;
                }
                ++this_end_index;
            }
        }
        globals.cx1.op_[t].op_start_index = this_start_index;
        globals.cx1.op_[t].op_end_index = this_end_index;
        last_end_index = this_end_index;
    }

    // last partition
    globals.cx1.op_[globals.num_output_threads - 1].op_start_index = last_end_index;
    globals.cx1.op_[globals.num_output_threads - 1].op_end_index = globals.lv2_num_items_db;

    memset(globals.lv2_aux, 0, sizeof(globals.lv2_aux[0]) * globals.lv2_num_items_db);
    pthread_barrier_init(&globals.output_barrier, NULL, globals.num_output_threads);
}

void* lv2_output(void* _op) {
    outputpartition_data_t *op = (outputpartition_data_t*) _op;
    edge2sdbg_global_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int start_idx, end_idx;
    int has_solid_a = 0; // has solid (k+1)-mer aSb
    int has_solid_b = 0; // has solid aSb
    int last_a[4], outputed_b;

    for (start_idx = op_start_index; start_idx < op_end_index; start_idx = end_idx) {
        end_idx = start_idx + 1;
        edge_word_t *item = globals.lv2_substrings_db + globals.permutation_db[start_idx];
        while (end_idx < op_end_index && 
               !IsDiffKMinusOneMer(
                    item, 
                    globals.lv2_substrings_db + globals.permutation_db[end_idx],
                    globals.lv2_num_items_db,
                    globals.kmer_k)) {
            ++end_idx;
        }

        // clean marking
        has_solid_a = has_solid_b = 0;
        outputed_b = 0;
        for (int i = start_idx; i < end_idx; ++i) {
            edge_word_t *cur_item = globals.lv2_substrings_db + globals.permutation_db[i];
            int a = Extract_a(cur_item, globals.words_per_substring, globals.lv2_num_items_db, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, globals.lv2_num_items_db);

            if (a != kSentinelValue && b != kSentinelValue) {
                has_solid_a |= 1 << a;
                has_solid_b |= 1 << b;
            }
            if (a != kSentinelValue && 
                (b != kSentinelValue || !(has_solid_a & (1 << a)))) {
                last_a[a] = i;
            }
        }

        for (int i = start_idx, j; i < end_idx; i = j) {
            edge_word_t *cur_item = globals.lv2_substrings_db + globals.permutation_db[i];
            int a = Extract_a(cur_item, globals.words_per_substring, globals.lv2_num_items_db, globals.kmer_k);
            int b = Extract_b(cur_item, globals.words_per_substring, globals.lv2_num_items_db);

            j = i + 1;
            while (j < end_idx) {
                edge_word_t *next_item = globals.lv2_substrings_db + globals.permutation_db[j];
                if (Extract_a(next_item, globals.words_per_substring, globals.lv2_num_items_db, globals.kmer_k) != a ||
                    Extract_b(next_item, globals.words_per_substring, globals.lv2_num_items_db) != b) {
                    break;
                } else {
                    ++j;
                }
            }

            int w, last, is_dollar = 0;

            if (a == kSentinelValue) {
                assert(b != kSentinelValue);
                if (has_solid_b & (1 << b)) {
                    continue;
                }
                is_dollar = 1;
            }

            if (b == kSentinelValue) {
                assert(a != kSentinelValue);
                if (has_solid_a & (1 << a)) {
                    continue;
                }
            }

            w = (b == kSentinelValue) ? 0 : ((outputed_b & (1 << b)) ? b + 5 : b + 1);
            outputed_b |= 1 << b;
            last = (a == kSentinelValue) ? 0 : ((last_a[a] == j - 1) ? 1 : 0);

            assert(!(globals.lv2_aux[i] & (1 << 7)));
            globals.lv2_aux[i] = w | (last << 4) | (is_dollar << 5) | (1 << 7);
        }
    }

    pthread_barrier_wait(&globals.output_barrier);

    if (op_start_index == 0) {
        xtimer_t local_timer;
        local_timer.reset();
        local_timer.start();
        for (int i = 0; i < globals.lv2_num_items_db; ++i) {
            if (globals.lv2_aux[i] & (1 << 7)) {
                edge_word_t *item = globals.lv2_substrings_db + globals.permutation_db[i];
                while (ExtractFirstChar(item) > globals.cur_suffix_first_char) {
                    ++globals.cur_suffix_first_char;
                    fprintf(globals.output_f_file, "%lld\n", (long long)globals.total_number_edges);
                }

                multi_t counting_db = std::min(kMaxMulti_t, 
                    ExtractCounting(item, globals.words_per_substring, globals.lv2_num_items_db));
                // output
                globals.sdbg_writer.outputW(globals.lv2_aux[i] & 0xF);
                globals.sdbg_writer.outputLast((globals.lv2_aux[i] >> 4) & 1);
                globals.sdbg_writer.outputIsDollar((globals.lv2_aux[i] >> 5) & 1);
                if (counting_db <= kMaxMulti2_t) {
                    multi2_t c = counting_db;
                    fwrite(&c, sizeof(multi2_t), 1, globals.output_multiplicity_file);   
                } else {
                    int64_t c = counting_db | (globals.total_number_edges << 16);
                    fwrite(&c, sizeof(int64_t), 1, globals.output_multiplicity_file2);
                    fwrite(&kMulti2Sp, sizeof(multi2_t), 1, globals.output_multiplicity_file);
                }

                globals.total_number_edges++;
                globals.num_chars_in_w[globals.lv2_aux[i] & 0xF]++;
                globals.num_ones_in_last += (globals.lv2_aux[i] >> 4) & 1;

                if ((globals.lv2_aux[i] >> 5) & 1) {
                    globals.num_dollar_nodes++;
                    if (globals.num_dollar_nodes >= kMaxDummyEdges) {
                        err("[ERROR B::%s] Too many dummy nodes (>= %lld)! The graph contains too many tips!\n", __func__, (long long)kMaxDummyEdges);
                        exit(1);
                    }
                    for (int64_t i = 0; i < globals.words_per_dummy_node; ++i) {
                        globals.dummy_nodes_writer.output(item[i * globals.lv2_num_items_db]);
                    }
                }
                if ((globals.lv2_aux[i] & 0xF) == 0) {
                    globals.num_dummy_edges++;
                }
            }
        }
        local_timer.stop();

        if (cx1_t::kCX1Verbose >= 4) {
            log("[B::%s] SdBG calc linear part: %lf\n", __func__, local_timer.elapsed());
        }
    }
    return NULL;
}

void lv2_post_output(edge2sdbg_global_t &globals) {
    pthread_barrier_destroy(&globals.output_barrier);
}

void post_proc(edge2sdbg_global_t &globals) {
	if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] Number of $ A C G T A- C- G- T-:\n", __func__);
    }

    for (int i = 0; i < 9; ++i) {
        log("%lld ", globals.num_chars_in_w[i]);
    }
    log("\n");

    // --- write tails ---
    fprintf(globals.output_f_file, "%lld\n", (long long)globals.total_number_edges);
    fprintf(globals.output_f_file, "%d\n", globals.kmer_k);
    fprintf(globals.output_f_file, "%lld\n", (long long)globals.num_dollar_nodes);

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] Total number of edges: %llu\n", __func__, globals.total_number_edges);
        log("[B::%s] Total number of ONEs: %llu\n", __func__, globals.num_ones_in_last);
        log("[B::%s] Total number of v$ edges: %llu\n", __func__, globals.num_dummy_edges);
        log("[B::%s] Total number of $v edges: %llu\n", __func__, globals.num_dollar_nodes);
    }

    // --- cleaning ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.packed_edges);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.lv2_substrings_db);
    free(globals.permutation_db);
    free(globals.lv2_aux);
    fclose(globals.output_f_file);
    fclose(globals.output_multiplicity_file);
    fclose(globals.output_multiplicity_file2);
    if (globals.mult_mem_type == 1) {
        free(globals.multiplicity8);
    } else if (globals.mult_mem_type == 2) {
        free(globals.multiplicity16);
    }
    globals.dummy_nodes_writer.destroy();
 #ifndef USE_GPU
    free(globals.cpu_sort_space);
 #else
    free_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);
 #endif
}

} // namespace