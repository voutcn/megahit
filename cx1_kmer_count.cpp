#include "cx1_kmer_count.h"

#include <string.h>
#include <algorithm>
#include <zlib.h>
#include <omp.h>

#include "sdbg_builder_writers.h"
#include "mem_file_checker-inl.h"
#include "kseq.h"
#include "utils.h"
#include "kmer_uint32.h"

#ifdef DISABLE_GPU
#include "lv2_cpu_sort.h"
#else
#include "lv2_gpu_functions.h"
#endif 

namespace cx1_kmer_count {

// helpers
KSEQ_INIT(gzFile, gzread)
typedef CX1<count_global_t, kNumBuckets> cx1_t;
typedef CX1<count_global_t, kNumBuckets>::readpartition_data_t readpartition_data_t;
typedef CX1<count_global_t, kNumBuckets>::bucketpartition_data_t bucketpartition_data_t;
typedef CX1<count_global_t, kNumBuckets>::outputpartition_data_t outputpartition_data_t;

int dna_map[256];

/**
 * @brief encode read_id and its offset in one int64_t
 */
inline int64_t EncodeOffset(int64_t read_id, int offset, int strand, int length_num_bits) {
    return (read_id << (length_num_bits + 1)) | (offset << 1) | strand;
}

/*
 * Packs an ASCII read into 2-bit per base form. The last base goes to the MSB of the first word.
 * -Params-
 * read: the read, in ASCII ACGT
 * p: a pointer to the first edge_word_t of the packed sequence to be written to
 * read_length: number of bases in the read
 * last_word_shift: the number of empty bits in the last word. we need to shift the bits up by this amount. solely determined by read_length.
 */
inline void PackReadFromAscii(char* read, edge_word_t* p, int read_length, int words_per_read) {
    // for de Bruijn graph construction, packing the reverse is more convenient
    edge_word_t w = 0;
    int i, j;
    for (i = 0, j = 0; j < read_length; ++j) {
        if (j % kCharsPerEdgeWord == 0 && j) { // TODO bitwise?
            *(p++) = w;
            w = 0;
            ++i;
        }
        while (read[j] == 'N') {
            break;
        }
        w = (w << kBitsPerEdgeChar) | dna_map[ (int)read[ j ] ];
    }

    int last_word_shift = j % kCharsPerEdgeWord;
    last_word_shift = last_word_shift ? (kCharsPerEdgeWord - last_word_shift) * kBitsPerEdgeChar : 0;
    *p = w << last_word_shift;

    while (++i < words_per_read) {
        *(++p) = 0;
    }

    *p |= read_length;
}

inline edge_word_t* GetReadPtr(int64_t i, count_global_t &globals) {
	return globals.packed_reads + i * globals.words_per_read;
}

inline int GetReadLength(edge_word_t* read_p, int words_per_read, int mask) {
    return *(read_p + words_per_read - 1) & mask;
}

inline int GetReadLengthByID(int64_t id, count_global_t &globals) {
    return *(globals.packed_reads + (id + 1) * globals.words_per_read - 1) & globals.read_length_mask;
}

/**
 * @brief extract the nth char in a packed read/edge
 */
inline int ExtractNthChar(edge_word_t *read_ptr, int n) {
    int which_word = n / kCharsPerEdgeWord;
    int index_in_word = n % kCharsPerEdgeWord;
    return (read_ptr[which_word] >> (kBitsPerEdgeChar * (kCharsPerEdgeWord - 1 - index_in_word))) & kEdgeCharMask;
}

// 'spacing' is the strip length for read-word "coalescing"
void CopySubstring(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, count_global_t &globals) {
    int64_t spacing = globals.cx1.lv2_num_items_;
    int words_per_read = globals.words_per_read;
    int words_per_substring = globals.words_per_substring;

    // copy words of the suffix to the suffix pool
    int which_word = offset / kCharsPerEdgeWord;
    int word_offset = offset % kCharsPerEdgeWord;
    edge_word_t *src_p = src_read + which_word;
    edge_word_t *dest_p = dest;
    int num_words_copied = 0;
    if (!word_offset) { // special case (word aligned), easy
        while (which_word < words_per_read && num_words_copied < words_per_substring) {
            *dest_p = *src_p; // write out
            dest_p += spacing;
            src_p++;
            which_word++;
            num_words_copied++;
        }
    } else { // not word-aligned
        int bit_shift = offset * kBitsPerEdgeChar;
        edge_word_t s = *src_p;
        edge_word_t d = s << bit_shift;
        which_word++;
        while (which_word < words_per_read) {
            s = *(++src_p);
            d |= s >> (kBitsPerEdgeWord - bit_shift);
            *dest_p = d; // write out
            if (++num_words_copied >= words_per_substring) goto here;
            dest_p += spacing;
            d = s << bit_shift;
            which_word++;
        }
        *dest_p = d; // write last word
	here:
        ;
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < globals.words_per_substring) {
            *p = 0;
        }
        which_word++;
        while (which_word < globals.words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }
}

void CopySubstringRC(edge_word_t* dest, edge_word_t* src_read, int offset, int num_chars_to_copy, count_global_t &globals) {
    assert(num_chars_to_copy == globals.kmer_k + 1);
    int spacing = globals.cx1.lv2_num_items_;
    int which_word = (offset + num_chars_to_copy - 1) / kCharsPerEdgeWord;
    int word_offset = (offset + num_chars_to_copy - 1) % kCharsPerEdgeWord;
    edge_word_t *dest_p = dest;

    if (word_offset == kCharsPerEdgeWord - 1) { // edge_word_t aligned
        for (int i = 0; i < globals.words_per_substring && i <= which_word; ++i) {
            *dest_p = ~ mirror(src_read[which_word - i]);
            dest_p += spacing;
        }
    } else {
        int bit_offset = (kCharsPerEdgeWord - 1 - word_offset) * kBitsPerEdgeChar;
        int i;
        edge_word_t w;
        for (i = 0; i < globals.words_per_substring - 1 && i < which_word; ++i) {
            w = (src_read[which_word - i] >> bit_offset) |
                                      (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
            *dest_p = ~ mirror(w);
            dest_p += spacing;
        }
        // last word
        w = src_read[which_word - i] >> bit_offset;
        if (which_word >= i + 1) {
            w |= (src_read[which_word - i - 1] << (kBitsPerEdgeWord - bit_offset));
        }
        *dest_p = ~ mirror(w);
    }

    {
        // now mask the extra bits (TODO can be optimized)
        int num_bits_to_copy = num_chars_to_copy * 2;
        int which_word = num_bits_to_copy / kBitsPerEdgeWord;
        edge_word_t *p = dest + which_word * spacing;
        int bits_to_clear = kBitsPerEdgeWord - num_bits_to_copy % kBitsPerEdgeWord;
        if (bits_to_clear < kBitsPerEdgeWord) {
            *p >>= bits_to_clear;
            *p <<= bits_to_clear;
        } else if (which_word < globals.words_per_substring) {
            *p = 0;
        }
        which_word++;
        while (which_word < globals.words_per_substring) { // fill zero
            *(p+=spacing) = 0;
            which_word++;
        }
    }
}

inline bool IsDifferentEdges(edge_word_t *item1, edge_word_t* item2, int num_words, int spacing) {
    for (int i = num_words - 1; i >= 0; --i) {
        if (*(item1 + (int64_t)i * spacing) != *(item2 + (int64_t)i * spacing)) {
            return true;
        }
    }
    return false;
}

/**
 * @brief pack an edge and its multiplicity to word-aligned spaces
 */
inline void PackEdge(edge_word_t *dest, edge_word_t *item, int counting, struct count_global_t &globals) {
    for (int i = 0; i < globals.words_per_edge && i < globals.words_per_substring; ++i) {
        dest[i] = *(item + (int64_t)i * globals.lv2_num_items_db);
    }
    int chars_in_last_word = (globals.kmer_k + 1) % kCharsPerEdgeWord;
    int which_word = (globals.kmer_k + 1) / kCharsPerEdgeWord;
    if (chars_in_last_word > 0) {
        dest[which_word] >>= (kCharsPerEdgeWord - chars_in_last_word) * 2;
        dest[which_word] <<= (kCharsPerEdgeWord - chars_in_last_word) * 2;
    } else {
        dest[which_word] = 0;
    }
    while (++which_word < globals.words_per_edge) {
        dest[which_word] = 0;
    }

    dest[globals.words_per_edge - 1] |= std::min(kMaxMulti_t, counting);
}

// function pass to CX1

int64_t encode_lv1_diff_base(int64_t read_id, count_global_t &g) {
	return EncodeOffset(read_id, 0, 0, g.offset_num_bits);
}

void read_input_prepare(count_global_t &globals) { // num_items_, num_cpu_threads_ and num_output_threads_ must be set here
	for (int i = 0; i < 10; ++i) {
		dna_map[int("ACGTNacgtn"[i])] = "0123101231"[i] - '0';
	}

	int64_t num_reads = 0;
    int bits_read_length = 1; // bit needed to store read_length
    while ((1 << bits_read_length) - 1 < globals.max_read_length) {
        ++bits_read_length;
    }
    globals.read_length_mask = (1 << bits_read_length) - 1;

    {
        globals.offset_num_bits = 0;
        int len = 1;
        while (len - 1 < globals.max_read_length) {
            globals.offset_num_bits++;
            len *= 2;
        }
    }

    globals.words_per_read = DivCeiling(globals.max_read_length * kBitsPerEdgeChar + bits_read_length, kBitsPerEdgeWord);
    int64_t max_num_reads = globals.host_mem / (sizeof(edge_word_t) * globals.words_per_read) * 3 / 4; //TODO: more accurate
    int read_length;
    edge_word_t *packed_reads;
    edge_word_t *packed_reads_p; // current pointer
    globals.capacity = std::min(max_num_reads, int64_t(1048576)); // initial capacity 1M
    gzFile fp = strcmp(globals.input_file, "-") ? gzopen(globals.input_file, "r") : gzdopen(fileno(stdin), "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    packed_reads_p = packed_reads = (edge_word_t*) MallocAndCheck(globals.capacity * globals.words_per_read * sizeof(edge_word_t), __FILE__, __LINE__);
    if (cx1_t::kCX1Verbose >= 2) {
        log("[C::%s] Max read length is %d; words per read: %d\n", __func__, globals.max_read_length, globals.words_per_read);
    }

    // --- main reading loop ---
    bool stop_reading = false;
    while ((read_length = kseq_read(seq)) >= 0 && !stop_reading) {
        std::reverse(seq->seq.s, seq->seq.s + read_length);
        char *next_p = seq->seq.s;
        while (read_length > globals.kmer_k) {
            int scan_len = 0;
            while (scan_len < read_length && next_p[scan_len] != 'N') {
                ++scan_len;
            }

            if (scan_len > globals.kmer_k && scan_len <= globals.max_read_length) {
                if (num_reads >= globals.capacity) {
                    if (globals.capacity == max_num_reads) {
                        err("[C::%s WRANING] No enough memory to hold all the reads. Only the first %llu reads are kept.\n", __func__, num_reads);
                        stop_reading = true;
                        break;
                    } 
                    globals.capacity = std::min(globals.capacity * 2, max_num_reads);
                    edge_word_t *new_ptr = (edge_word_t*) realloc(packed_reads, globals.capacity * globals.words_per_read * sizeof(edge_word_t));
                    if (new_ptr != NULL) {
                        packed_reads = new_ptr;
                        packed_reads_p = packed_reads + globals.words_per_read * num_reads;
                        globals.capacity = globals.capacity;
                    } else {
                        err("[C::%s WRANING] No enough memory to hold all the reads. Only the first %llu reads are kept.\n", __func__, num_reads);
                        stop_reading = true;
                        break;
                    }
                }
                // read length is ok! compress and store the packed read
                PackReadFromAscii(next_p, packed_reads_p, scan_len, globals.words_per_read);
                packed_reads_p += globals.words_per_read;
                ++num_reads;
            } else if (scan_len > globals.max_read_length) { // this read length is wrong
                err("[C::%s WARNING] Found a read of length %d > max read length = %d\n, it will be discarded.", __func__, scan_len, globals.max_read_length);
            }

            while (scan_len < read_length && next_p[scan_len] == 'N') {
                ++scan_len;
            }
            read_length -= scan_len;
            next_p += scan_len;
        }
    }

    globals.num_reads = num_reads;
    globals.mem_packed_reads = globals.num_reads * globals.words_per_read * sizeof(edge_word_t);
    globals.packed_reads = (edge_word_t*) ReAllocAndCheck(packed_reads, globals.mem_packed_reads, __FILE__, __LINE__);
    if (!globals.packed_reads) {
        err("[C::%s ERROR] Cannot reallocate memory for packed reads!\n", __func__);
        exit(1);
    }
    log("[C::%s] Total number of reads: %lld\n", __func__, (long long)globals.num_reads);

    // set cx1 param
    globals.cx1.num_cpu_threads_ = globals.num_cpu_threads;
    globals.cx1.num_output_threads_ = globals.num_output_threads;
    globals.cx1.num_items_ = globals.num_reads;

    kseq_destroy(seq);
    gzclose(fp);
}

void* lv0_calc_bucket_size(void* _data) {
	readpartition_data_t &rp = *((readpartition_data_t*) _data);
    count_global_t &globals = *(rp.globals);
    int64_t *bucket_sizes = rp.rp_bucket_sizes;
    memset(bucket_sizes, 0, kNumBuckets * sizeof(int64_t));
    edge_word_t *read_p = GetReadPtr(rp.rp_start_id, globals);
    KmerUint32 edge, rev_edge; // (k+1)-mer and its rc
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, read_p += globals.words_per_read) {
        int read_length = GetReadLength(read_p, globals.words_per_read, globals.read_length_mask);
        if (read_length < globals.kmer_k + 1) { continue; }
        edge.init(read_p, globals.kmer_k + 1);
        rev_edge.clean();
        for (int i = 0; i <= globals.kmer_k; ++i) {
            rev_edge.Append(3 - ExtractNthChar(read_p, globals.kmer_k - i));
        }

        int last_char_offset = globals.kmer_k;
        while (true) {
            if (rev_edge < edge) {
                bucket_sizes[rev_edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            } else {
                bucket_sizes[edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar]++;
            }

            if (++last_char_offset >= read_length) {
                break;
            } else {
                int c = ExtractNthChar(read_p, last_char_offset);
                edge.ShiftLeftAppend(c);
                rev_edge.ShiftRightAppend(3 - c);
            }
        }
    }
    return NULL;
}

void init_global_and_set_cx1(count_global_t &globals) {
    // --- calculate lv2 memory ---
    globals.max_bucket_size = *std::max_element(globals.cx1.bucket_sizes_, globals.cx1.bucket_sizes_ + kNumBuckets);
    globals.tot_bucket_size = 0;
    for (int i = 0; i < kNumBuckets; ++i) { globals.tot_bucket_size += globals.cx1.bucket_sizes_[i]; }

#ifdef DISABLE_GPU
    globals.cx1.max_lv2_items_ = std::max(globals.max_bucket_size, kMinLv2BatchSize);
#else
    int64_t lv2_mem = globals.gpu_mem - 1073741824; // should reserver ~1G for GPU sorting
    globals.cx1.max_lv2_items_ = lv2_mem / cx1_t::kGPUBytePerItem;
    globals.cx1.max_lv2_items_ = std::min(globals.cx1.max_lv2_items_, globals.tot_bucket_size);
    if (globals.max_bucket_size > globals.cx1.max_lv2_items_) {
        err("[ERROR C::%s] Bucket too large for GPU: contains %lld items. Please try CPU version.\n", __func__, globals.max_bucket_size);
        // TODO: auto switch to CPU version
        exit(1);
    }
    globals.cx1.max_lv2_items_ = std::max(globals.max_bucket_size, kMinLv2BatchSizeGPU);
#endif
    globals.words_per_substring = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar, kBitsPerEdgeWord);
    globals.words_per_edge = DivCeiling((globals.kmer_k + 1) * kBitsPerEdgeChar + kBitsPerMulti_t, kBitsPerEdgeWord);
    // lv2 bytes: substring, permutation, readinfo
    int64_t lv2_bytes_per_item = (globals.words_per_substring) * sizeof(edge_word_t) + sizeof(uint32_t) + sizeof(int64_t);
    lv2_bytes_per_item = lv2_bytes_per_item * 2; // double buffering
#ifdef DISABLE_GPU
    lv2_bytes_per_item += sizeof(uint64_t) * 2; // CPU memory is used to simulate GPU
#endif

    if (cx1_t::kCX1Verbose >= 2) {
        log("[C::%s] %d words per substring, %d words per edge\n", __func__, globals.words_per_substring, globals.words_per_edge);
    }
    // --- memory stuff ---
    int64_t mem_remained = globals.host_mem 
                         - globals.mem_packed_reads
                         - globals.num_reads * sizeof(unsigned char) * 2 // first_in0 & last_out0
                         - kNumBuckets * sizeof(int64_t) * (globals.num_cpu_threads * 3 + 1)
                         - (kMaxMulti_t + 1) * (globals.num_output_threads + 1) * sizeof(int64_t);
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
        log("[C::%s] Memory for reads: %lld\n", __func__, globals.mem_packed_reads);
        log("[C::%s] max # lv.1 items = %lld\n", __func__, globals.cx1.max_lv1_items_);
        log("[C::%s] max # lv.2 items = %lld\n", __func__, globals.cx1.max_lv2_items_);
    }

    // --- alloc memory ---
    globals.lv1_items = (int*) MallocAndCheck(globals.cx1.max_lv1_items_ * sizeof(int), __FILE__, __LINE__);
    globals.lv2_substrings = (edge_word_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.lv2_substrings_db = (edge_word_t*) MallocAndCheck(globals.cx1.max_lv2_items_ * globals.words_per_substring * sizeof(edge_word_t), __FILE__, __LINE__);
    globals.permutation = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.permutation_db = (uint32_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(uint32_t), __FILE__, __LINE__);
    globals.lv2_read_info = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
    globals.lv2_read_info_db = (int64_t *) MallocAndCheck(globals.cx1.max_lv2_items_ * sizeof(int64_t), __FILE__, __LINE__);
#ifdef LONG_READS
    globals.first_0_out = (uint16_t*) MallocAndCheck(globals.num_reads * sizeof(uint16_t), __FILE__, __LINE__);
    globals.last_0_in = (uint16_t*) MallocAndCheck(globals.num_reads * sizeof(uint16_t), __FILE__, __LINE__);
#else
    globals.first_0_out = (unsigned char*) MallocAndCheck(globals.num_reads * sizeof(unsigned char), __FILE__, __LINE__);
    globals.last_0_in = (unsigned char*) MallocAndCheck(globals.num_reads * sizeof(unsigned char), __FILE__, __LINE__);
#endif
 #ifdef DISABLE_GPU
    globals.cpu_sort_space = (uint64_t*) MallocAndCheck(sizeof(uint64_t) * globals.cx1.max_lv2_items_, __FILE__, __LINE__);
 #else
    alloc_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2, (size_t)globals.cx1.max_lv2_items_);
 #endif
    memset(globals.first_0_out, 0xFF, globals.num_reads * sizeof(globals.first_0_out[0]));
    memset(globals.last_0_in, 0xFF, globals.num_reads * sizeof(globals.last_0_in[0]));

    // --- initialize stat ---
    globals.edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * sizeof(int64_t), __FILE__, __LINE__);
    globals.thread_edge_counting = (int64_t *) MallocAndCheck((kMaxMulti_t + 1) * globals.num_output_threads * sizeof(int64_t), __FILE__, __LINE__);
    memset(globals.edge_counting, 0, (kMaxMulti_t + 1) * sizeof(int64_t));

    // --- initialize lock ---
    pthread_mutex_init(&globals.lv1_items_scanning_lock, NULL); // init lock

    // --- initialize writer ---
    globals.word_writer = (WordWriter*) MallocAndCheck(globals.num_output_threads * sizeof(WordWriter), __FILE__, __LINE__);
    for (int t = 0; t < globals.num_output_threads; ++t) {
        char edges_file_name[10240];
        sprintf(edges_file_name, "%s.edges.%d", globals.output_prefix, t);
        globals.word_writer[t].init(edges_file_name);
    }

    // --- write the edge file header ---
    globals.word_writer[0].output(globals.kmer_k);
    globals.word_writer[0].output(globals.words_per_edge);
}

void* lv1_fill_offset(void* _data) {
	readpartition_data_t &rp = *((readpartition_data_t*) _data);
    count_global_t &globals = *(rp.globals);
    int64_t *prev_full_offsets = (int64_t *)MallocAndCheck(kNumBuckets * sizeof(int64_t), __FILE__, __LINE__); // temporary array for computing differentials
    assert(prev_full_offsets != NULL);
    for (int b = globals.cx1.lv1_start_bucket_; b < globals.cx1.lv1_end_bucket_; ++b)
        prev_full_offsets[b] = rp.rp_lv1_differential_base;
    // this loop is VERY similar to that in PreprocessScanToFillBucketSizesThread
    edge_word_t *read_p = GetReadPtr(rp.rp_start_id, globals);
    KmerUint32 edge, rev_edge; // (k+1)-mer and its rc
    int key;
    for (int64_t read_id = rp.rp_start_id; read_id < rp.rp_end_id; ++read_id, read_p += globals.words_per_read) {
        int read_length = GetReadLength(read_p, globals.words_per_read, globals.read_length_mask);
        if (read_length < globals.kmer_k + 1) { continue; }
        edge.init(read_p, globals.kmer_k + 1);
        rev_edge.clean();
        for (int i = 0; i <= globals.kmer_k; ++i) {
            rev_edge.Append(3 - ExtractNthChar(read_p, globals.kmer_k - i));
        }

        // ===== this is a macro to save some copy&paste ================
#define CHECK_AND_SAVE_OFFSET(offset, strand)    \
    do {                                                                \
      assert(offset + globals.kmer_k < read_length); \
      if (((key - globals.cx1.lv1_start_bucket_) ^ (key - globals.cx1.lv1_end_bucket_)) & kSignBitMask) { \
        int64_t full_offset = EncodeOffset(read_id, offset, strand, globals.offset_num_bits); \
        int64_t differential = full_offset - prev_full_offsets[key];      \
        if (differential > cx1_t::kDifferentialLimit) {                      \
          pthread_mutex_lock(&globals.lv1_items_scanning_lock); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = -globals.cx1.lv1_items_special_.size() - 1; \
          globals.cx1.lv1_items_special_.push_back(full_offset);                  \
          pthread_mutex_unlock(&globals.lv1_items_scanning_lock); \
        } else {                                                              \
          assert((int) differential >= 0); \
          globals.lv1_items[ rp.rp_bucket_offsets[key]++ ] = (int) differential; \
        } \
        prev_full_offsets[key] = full_offset;                           \
      }                                                                 \
    } while (0)
        // ^^^^^ why is the macro surrounded by a do-while? please ask Google
        // =========== end macro ==========================

        // shift the key char by char
        int last_char_offset = globals.kmer_k;
        while (true) {
            if (rev_edge < edge) {
                key = rev_edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 1);
            } else {
                key = edge.data_[0] >> (kCharsPerEdgeWord - kBucketPrefixLength) * kBitsPerEdgeChar;
                CHECK_AND_SAVE_OFFSET(last_char_offset - globals.kmer_k, 0);
            }

            if (++last_char_offset >= read_length) {
                break;
            } else {
                int c = ExtractNthChar(read_p, last_char_offset);
                edge.ShiftLeftAppend(c);
                rev_edge.ShiftRightAppend(3 - c);
            }
        }
    }

#undef CHECK_AND_SAVE_OFFSET

    free(prev_full_offsets);
    return NULL;
}

void* lv2_extract_substr(void* _data) {
	bucketpartition_data_t &bp = *((bucketpartition_data_t*) _data);
    count_global_t &globals = *(bp.globals);
    int *lv1_p = globals.lv1_items + globals.cx1.rp_[0].rp_bucket_offsets[ bp.bp_start_bucket ];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1; // 0000....00011..11
    edge_word_t *substrings_p = globals.lv2_substrings +
                         (globals.cx1.rp_[0].rp_bucket_offsets[ bp.bp_start_bucket ] - globals.cx1.rp_[0].rp_bucket_offsets[ globals.cx1.lv2_start_bucket_ ]);
    int64_t *read_info_p = globals.lv2_read_info + 
                       (globals.cx1.rp_[0].rp_bucket_offsets[ bp.bp_start_bucket ] - globals.cx1.rp_[0].rp_bucket_offsets[ globals.cx1.lv2_start_bucket_ ]);
    for (int b = bp.bp_start_bucket; b < bp.bp_end_bucket; ++b) {
        for (int t = 0; t < globals.num_cpu_threads; ++t) {
            int64_t full_offset = globals.cx1.rp_[t].rp_lv1_differential_base;
            int num = globals.cx1.rp_[t].rp_bucket_sizes[b];
            for (int i = 0; i < num; ++i) {
                if (*lv1_p >= 0) {
                    full_offset += *(lv1_p++);
                } else {
                    full_offset = globals.cx1.lv1_items_special_[-1 - *(lv1_p++)];
                }
                int64_t read_id = full_offset >> (globals.offset_num_bits + 1);
                int strand = full_offset & 1;
                int offset = (full_offset >> 1) & offset_mask;
                int num_chars_to_copy = globals.kmer_k + 1;
                unsigned char prev, next;
                if (offset > 0) {
                    prev = ExtractNthChar(GetReadPtr(read_id, globals), offset - 1);
                } else {
                    prev = kSentinelValue;
                }

                if (offset + globals.kmer_k + 1 < GetReadLengthByID(read_id, globals)) {
                    next = ExtractNthChar(GetReadPtr(read_id, globals), offset + globals.kmer_k + 1);
                } else {
                    next = kSentinelValue;
                }

                if (strand == 0) {
                    CopySubstring(substrings_p, GetReadPtr(read_id, globals), offset, num_chars_to_copy, globals);
                    *read_info_p = (full_offset << 6) | (prev << 3) | next;
                } else {
                    CopySubstringRC(substrings_p, GetReadPtr(read_id, globals), offset, num_chars_to_copy, globals);
                    *read_info_p = (full_offset << 6) | ((next == kSentinelValue ? kSentinelValue : (3 - next)) << 3)
                                                      | (prev == kSentinelValue ? kSentinelValue : (3 - prev));
                }
                substrings_p++;
                read_info_p++;
            }
        }
    }
    return NULL;
}

void lv2_sort(count_global_t &globals) {
	xtimer_t local_timer;
#ifdef DISABLE_GPU
    if (cx1_t::kCX1Verbose >= 4) {
	    local_timer.reset();
	    local_timer.start();
    }
    omp_set_num_threads(globals.num_cpu_threads - globals.num_output_threads);
    lv2_cpu_sort(globals.lv2_substrings, globals.permutation, globals.cpu_sort_space, globals.words_per_substring, globals.cx1.lv2_num_items_);
    omp_set_num_threads(globals.num_cpu_threads);
    local_timer.stop();

    if (cx1_t::kCX1Verbose >= 4) {
        log("[C::%s] Sorting substrings with CPU...done. Time elapsed: %.4lf\n", __func__, local_timer.elapsed());
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
        log("[C::%s] Sorting substrings with GPU...done. Time elapsed: %.4lf\n", __func__, local_timer.elapsed());
    }
#endif
}

void lv2_pre_output_partition(count_global_t &globals) {
	// swap double buffers
    globals.lv2_num_items_db = globals.cx1.lv2_num_items_;
    std::swap(globals.lv2_substrings_db, globals.lv2_substrings);
    std::swap(globals.permutation_db, globals.permutation);
    std::swap(globals.lv2_read_info_db, globals.lv2_read_info);

    // distribute partition
    int64_t last_end_index = 0;
    int64_t items_per_thread = globals.lv2_num_items_db / globals.num_output_threads;

    for (int t = 0; t < globals.num_output_threads - 1; ++t) {
        int64_t this_start_index = last_end_index;
        int64_t this_end_index = this_start_index + items_per_thread;
        if (this_end_index > globals.lv2_num_items_db) { this_end_index = globals.lv2_num_items_db; }
        if (this_end_index > 0) {
            while (this_end_index < globals.lv2_num_items_db) {
                edge_word_t *prev_item = globals.lv2_substrings_db + (globals.permutation_db[this_end_index - 1]);
                edge_word_t *item = globals.lv2_substrings_db + (globals.permutation_db[this_end_index]);
                if (IsDifferentEdges(prev_item, item, globals.words_per_substring, globals.lv2_num_items_db)) {
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

    memset(globals.thread_edge_counting, 0, sizeof(int64_t) * (kMaxMulti_t + 1) * globals.num_output_threads);
}

void* lv2_output(void* _op) {
    xtimer_t local_timer;
	if (cx1_t::kCX1Verbose >= 4) {
	    local_timer.start();
	    local_timer.reset();
    }

	outputpartition_data_t *op = (outputpartition_data_t*) _op;
    count_global_t &globals = *(op->globals);
    int64_t op_start_index = op->op_start_index;
    int64_t op_end_index = op->op_end_index;
    int thread_id = op->op_id;
    int start_idx;
    int end_idx;
    edge_word_t packed_edge[32];
    int count_prev[5], count_next[5];
    int64_t offset_mask = (1 << globals.offset_num_bits) - 1;
    int64_t *thread_edge_counting = globals.thread_edge_counting + thread_id * (kMaxMulti_t + 1);

    for (int i = op_start_index; i < op_end_index; i = end_idx) {
        start_idx = i;
        end_idx = i + 1;
        edge_word_t *first_item = globals.lv2_substrings_db + (globals.permutation_db[i]);
        while (end_idx < op_end_index) {
            if (IsDifferentEdges(first_item,
                               globals.lv2_substrings_db + globals.permutation_db[end_idx],
                               globals.words_per_substring, globals.lv2_num_items_db)) {
                break;
            }
            ++end_idx;
        }
        int count = end_idx - start_idx;

        // update read's first and last

        memset(count_prev, 0, sizeof(int) * 4);
        memset(count_next, 0, sizeof(int) * 4);
        bool has_in = false;
        bool has_out = false;
        for (int j = start_idx; j < end_idx; ++j) {
            int prev_and_next = globals.lv2_read_info_db[globals.permutation_db[j]] & ((1 << 6) - 1);
            count_prev[prev_and_next >> 3]++;
            count_next[prev_and_next & 7]++;
        }

        for (int j = 0; j < 4; ++j) {
            if (count_prev[j] >= globals.kmer_freq_threshold) { has_in = true; }
            if (count_next[j] >= globals.kmer_freq_threshold) { has_out = true; }
        }

        if (!has_in && count >= globals.kmer_freq_threshold) {
            for (int j = start_idx; j < end_idx; ++j) {
                int64_t read_info = globals.lv2_read_info_db[globals.permutation_db[j]] >> 6;
                int strand = read_info & 1;
                int offset = (read_info >> 1) & offset_mask;
                int64_t read_id = read_info >> (1 + globals.offset_num_bits);

                if (strand == 0) {
                    // update last
                    while (true) {
                        auto old_value = globals.last_0_in[read_id];
                        if (old_value != kSentinelOffset && old_value >= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.last_0_in + read_id, old_value, offset)) {
                            break;
                        }
                    }
                } else {
                    // update first
                    offset++;
                    while (true) {
                        auto old_value = globals.first_0_out[read_id];
                        if (old_value <= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.first_0_out + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
            }
        }

        if (!has_out && count >= globals.kmer_freq_threshold) {
            for (int j = start_idx; j < end_idx; ++j) {
                int64_t read_info = globals.lv2_read_info_db[globals.permutation_db[j]] >> 6;
                int strand = read_info & 1;
                int offset = (read_info >> 1) & offset_mask;
                int64_t read_id = read_info >> (1 + globals.offset_num_bits);

                if (strand == 0) {
                    // update first
                    offset++;
                    while (true) {
                        auto old_value = globals.first_0_out[read_id];
                        if (old_value <= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.first_0_out + read_id, old_value, offset)) {
                            break;
                        }
                    }
                } else {
                    // update last
                    while (true) {
                        auto old_value = globals.last_0_in[read_id];
                        if (old_value != kSentinelOffset && old_value >= offset) { break; }
                        if (__sync_bool_compare_and_swap(globals.last_0_in + read_id, old_value, offset)) {
                            break;
                        }
                    }
                }
            }
        }

        ++thread_edge_counting[std::min(count, kMaxMulti_t)];
        if (count >= globals.kmer_freq_threshold) {
            PackEdge(packed_edge, first_item, count, globals);
            for (int x = 0; x < globals.words_per_edge; ++x) {
                globals.word_writer[thread_id].output(packed_edge[x]);
            }
        }
    }

    if (cx1_t::kCX1Verbose >= 4) {
    	local_timer.stop();
        log("[C::%s] Counting time elapsed: %.4lfs\n", __func__, local_timer.elapsed());
    }
    return NULL;
}

void lv2_post_output(count_global_t &globals) {
    for (int t = 0; t < globals.num_output_threads; ++t) {
		for (int i = 1; i <= kMaxMulti_t; ++i) {
	        globals.edge_counting[i] += globals.thread_edge_counting[t * (kMaxMulti_t + 1) + i];
	    }
	}
}

void post_proc(count_global_t &globals) {
	// --- output reads for mercy ---
    int64_t num_candidate_reads = 0;
    int64_t num_has_tips = 0;
    FILE *candidate_file = OpenFileAndCheck((std::string(globals.output_prefix) + ".cand").c_str(), "wb");
    for (int64_t i = 0; i < globals.num_reads; ++i) {
        auto first = globals.first_0_out[i];
        auto last = globals.last_0_in[i];
        if (first != kSentinelOffset && last != kSentinelOffset) {
            ++num_has_tips;
            if (last > first) {
                ++num_candidate_reads;
                fwrite(GetReadPtr(i, globals), sizeof(uint32_t), globals.words_per_read, candidate_file);   
            }
        }
    }
    fclose(candidate_file);

    if (cx1_t::kCX1Verbose >= 2) {
        log("[C::%s] Total number of candidate reads: %lld(%lld)\n", __func__, num_candidate_reads, num_has_tips);
    }

    // --- stat ---
    int64_t num_solid_edges = 0;
    for (int i = globals.kmer_freq_threshold; i <= kMaxMulti_t; ++i) {
        num_solid_edges += globals.edge_counting[i];
    }

    if (cx1_t::kCX1Verbose >= 2) {
        log("[B::%s] Total number of solid edges: %llu\n", __func__, num_solid_edges);
    }

    FILE *counting_file = OpenFileAndCheck((std::string(globals.output_prefix)+".counting").c_str(), "w");
    for (int64_t i = 1, acc = 0; i <= kMaxMulti_t; ++i) {
        acc += globals.edge_counting[i];
        fprintf(counting_file, "%lld %lld\n", (long long)i, (long long)acc);
    }
    fclose(counting_file);

    // --- cleaning ---
    pthread_mutex_destroy(&globals.lv1_items_scanning_lock);
    free(globals.packed_reads);
    free(globals.lv1_items);
    free(globals.lv2_substrings);
    free(globals.permutation);
    free(globals.permutation_db);
    free(globals.lv2_substrings_db);
    free(globals.lv2_read_info);
    free(globals.lv2_read_info_db);
    free(globals.first_0_out);
    free(globals.last_0_in);
    free(globals.edge_counting);
    free(globals.thread_edge_counting);
    for (int t = 0; t < globals.num_output_threads; ++t) {        
        globals.word_writer[t].destroy();
    }
    free(globals.word_writer);

 #ifdef DISABLE_GPU
    free(globals.cpu_sort_space);
 #else
    free_gpu_buffers(globals.gpu_key_buffer1, globals.gpu_key_buffer2, globals.gpu_value_buffer1, globals.gpu_value_buffer2);
 #endif
}

} // namespace::cx1_kmer_count