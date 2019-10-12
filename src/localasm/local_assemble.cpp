#include "local_assemble.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>
#include "idba/contig_graph.h"
#include "idba/hash_graph.h"
#include "idba/sequence.h"
#include "kmlib/kmbit.h"

#include "hash_mapper.h"
#include "mapping_result_collector.h"
#include "sequence/io/contig/contig_reader.h"
#include "sequence/io/contig/contig_writer.h"
#include "sequence/io/sequence_lib.h"
#include "utils/histgram.h"
#include "utils/utils.h"

namespace {

static const int kMaxLocalRange = 650;
using TInsertSize = std::pair<double, double>;

void LaunchIDBA(const std::deque<Sequence> &reads, const Sequence &contig_end,
                std::deque<Sequence> &out_contigs,
                std::deque<ContigInfo> &out_contig_infos, uint32_t mink,
                uint32_t maxk, uint32_t step, ContigGraph &contig_graph) {
  int local_range = contig_end.size();
  HashGraph hash_graph;
  out_contigs.clear();
  out_contig_infos.clear();

  uint32_t max_read_len = 0;

  for (auto &read : reads) {
    max_read_len = std::max(max_read_len, read.size());
  }

  for (uint32_t kmer_size = mink; kmer_size <= std::min(maxk, max_read_len);
       kmer_size += step) {
    hash_graph.clear();
    hash_graph.set_kmer_size(kmer_size);

    for (auto &read : reads) {
      if (read.size() < kmer_size) continue;

      const Sequence seq(read);
      hash_graph.InsertKmers(seq);
    }

    auto histgram = hash_graph.coverage_histgram();
    double mean =
        histgram.percentile(1 - 1.0 * local_range / hash_graph.num_vertices());
    double threshold = mean;

    hash_graph.InsertKmers(contig_end);

    for (const auto &out_contig : out_contigs)
      hash_graph.InsertUncountKmers(out_contig);

    hash_graph.Assemble(out_contigs, out_contig_infos);

    contig_graph.clear();
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

std::vector<TInsertSize> EstimateInsertSize(
    const HashMapper &mapper, const SequenceLibCollection &lib_collection) {
  std::vector<TInsertSize> insert_sizes(lib_collection.size());
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto lib = lib_collection.GetLib(lib_id);

    if (!lib.IsPaired()) {
      continue;
    }

    Histgram<int> insert_hist;
    const size_t min_hist_size_for_estimation = 1u << 18;
    size_t processed_reads = 0;

    while (insert_hist.size() < min_hist_size_for_estimation &&
           processed_reads < lib.seq_count()) {
      size_t start_read_id = processed_reads;
      processed_reads = std::min(min_hist_size_for_estimation + start_read_id,
                                 lib.seq_count());

#pragma omp parallel for
      for (size_t i = start_read_id; i < processed_reads; i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        auto rec1 = mapper.TryMap(seq1);
        auto rec2 = mapper.TryMap(seq2);
        if (rec1.valid && rec2.valid) {
          if (rec1.contig_id == rec2.contig_id && rec1.strand != rec2.strand) {
            int insert_size;

            if (rec1.strand == 0) {
              insert_size = rec2.contig_to + seq2.length() - rec2.query_to -
                            (rec1.contig_from - rec1.query_from);
            } else {
              insert_size = rec1.contig_to + seq1.length() - rec1.query_to -
                            (rec2.contig_from - rec2.query_from);
            }

            if (insert_size >= (int)seq1.length() &&
                insert_size >= (int)seq2.length()) {
              insert_hist.insert(insert_size);
            }
          }
        }
      }
    }

    insert_hist.Trim(0.01);
    insert_sizes[lib_id] = TInsertSize(insert_hist.mean(), insert_hist.sd());

    xinfo("Lib {}, insert size: {.2} sd: {.2}\n", lib_id,
          insert_sizes[lib_id].first, insert_sizes[lib_id].second);
  }

  return insert_sizes;
}

int32_t LocalRange(const SequenceLib &lib, const TInsertSize &insert_size) {
  int32_t local_range = lib.GetMaxLength() - 1;

  if (lib.IsPaired() && insert_size.first >= lib.GetMaxLength()) {
    local_range = std::min(2 * insert_size.first,
                           insert_size.first + 3 * insert_size.second);
  }

  if (local_range > kMaxLocalRange) {
    local_range = kMaxLocalRange;
  }

  return local_range;
}

int32_t GetMaxLocalRange(const SequenceLibCollection &lib_collection,
                         const std::vector<TInsertSize> &insert_sizes) {
  int32_t max_local_range = 0;
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto &lib = lib_collection.GetLib(lib_id);
    max_local_range =
        std::max(max_local_range, LocalRange(lib, insert_sizes[lib_id]));
  }
  return max_local_range;
}

void MapToContigs(const HashMapper &mapper,
                  const SequenceLibCollection &lib_collection,
                  const std::vector<TInsertSize> &insert_sizes,
                  MappingResultCollector *collector) {
  for (unsigned lib_id = 0; lib_id < lib_collection.size(); ++lib_id) {
    auto &lib = lib_collection.GetLib(lib_id);
    int32_t local_range = LocalRange(lib, insert_sizes[lib_id]);
    bool is_paired = lib.IsPaired();

    size_t num_added = 0, num_mapped = 0;

    if (is_paired) {
#pragma omp parallel for reduction(+ : num_added, num_mapped)
      for (size_t i = 0; i < lib.seq_count(); i += 2) {
        auto seq1 = lib.GetSequenceView(i);
        auto seq2 = lib.GetSequenceView(i + 1);
        auto rec1 = mapper.TryMap(seq1);
        auto rec2 = mapper.TryMap(seq2);

        if (rec1.valid) {
          ++num_mapped;
          auto contig_len = mapper.refseq().GetSeqView(rec1.contig_id).length();
          num_added += collector->AddSingle(rec1, contig_len, seq1.length(),
                                            local_range);
          num_added += collector->AddMate(rec1, rec2, contig_len, seq2.id(),
                                          local_range);
        }

        if (rec2.valid) {
          ++num_mapped;
          auto contig_len = mapper.refseq().GetSeqView(rec2.contig_id).length();
          num_added += collector->AddSingle(rec2, contig_len, seq2.length(),
                                            local_range);
          num_added += collector->AddMate(rec2, rec1, contig_len, seq1.id(),
                                          local_range);
        }
      }
    } else {
#pragma omp parallel reduction(+ : num_added, num_mapped)
      for (size_t i = 0; i < lib.seq_count(); ++i) {
        auto seq = lib.GetSequenceView(i);
        auto rec = mapper.TryMap(seq);

        if (rec.valid) {
          ++num_mapped;
          num_added += collector->AddSingle(
              rec, mapper.refseq().GetSeqView(rec.contig_id).length(),
              seq.length(), local_range);
        }
      }
    }

    xinfo(
        "Lib {}: total {} reads, aligned {}, added {} reads to local "
        "assembly\n",
        lib_id, lib.seq_count(), num_mapped, num_added);
  }
}

void AssembleAndOutput(const HashMapper &mapper, const SeqPackage &read_pkg,
                       MappingResultCollector &result_collector,
                       const std::string &output_file,
                       const int32_t local_range,
                       const LocalAsmOption &opt) {
  const size_t min_num_reads = read_pkg.max_length() > 0 ?
      local_range / read_pkg.max_length(): 1;
  xinfo("Minimum number of reads to do local assembly: {}\n", min_num_reads);

  Sequence seq, contig_end;
  ContigGraph contig_graph;
  std::deque<Sequence> reads;
  std::deque<Sequence> out_contigs;
  std::deque<ContigInfo> out_contig_infos;

  ContigWriter local_contig_writer(output_file);

#pragma omp parallel for private(contig_graph, seq, contig_end, reads, \
                                 out_contigs, out_contig_infos)        \
    schedule(dynamic)
  for (uint64_t cid = 0; cid < mapper.refseq().seq_count(); ++cid) {
    auto contig_view = mapper.refseq().GetSeqView(cid);
    int cl = contig_view.length();

    for (uint8_t strand = 0; strand < 2; ++strand) {
      auto &mapping_rslts = result_collector.GetMappingResults(cid, strand);
      if (mapping_rslts.size() <= min_num_reads) {
        continue;
      }

      // collect local reads, convert them into Sequence
      reads.clear();
      uint64_t last_mapping_pos = -1;
      int pos_count = 0;

      for (const auto &encoded_rslt : mapping_rslts) {
        uint64_t pos = MappingResultCollector::GetContigAbsPos(encoded_rslt);
        pos_count = pos == last_mapping_pos ? pos_count + 1 : 1;
        last_mapping_pos = pos;

        if (pos_count <= 3) {
          seq.clear();
          auto read_view = read_pkg.GetSeqView(
              MappingResultCollector::GetReadId(encoded_rslt));

          for (unsigned ri = 0, rsz = read_view.length(); ri < rsz; ++ri) {
            seq.Append(read_view.base_at(ri));
          }
          reads.push_back(seq);
        }
      }

      contig_end.clear();

      if (strand == 0) {
        for (int j = 0, e = std::min(local_range, cl); j < e; ++j) {
          contig_end.Append(contig_view.base_at(j));
        }
      } else {
        for (int j = std::max(0, cl - local_range); j < cl; ++j) {
          contig_end.Append(contig_view.base_at(j));
        }
      }

      out_contigs.clear();
      LaunchIDBA(reads, contig_end, out_contigs, out_contig_infos, opt.kmin,
                 opt.kmax, opt.step, contig_graph);

      for (uint64_t j = 0; j < out_contigs.size(); ++j) {
        if (out_contigs[j].size() > opt.min_contig_len &&
            out_contigs[j].size() > opt.kmax) {
          auto str = out_contigs[j].str();
          local_contig_writer.WriteLocalContig(str, cid, strand, j);
        }
      }
    }
  }
}

}  // namespace

void RunLocalAssembly(const LocalAsmOption &opt) {
  SimpleTimer timer;
  timer.reset();
  timer.start();
  HashMapper mapper;
  mapper.LoadAndBuild(opt.contig_file, opt.min_contig_len, opt.seed_kmer,
                      opt.sparsity);
  mapper.SetMappingThreshold(opt.min_mapping_len, opt.similarity);
  timer.stop();
  xinfo("Hash mapper construction time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  SequenceLibCollection lib_collection;
  SeqPackage read_pkg;
  lib_collection.SetPath(opt.lib_file_prefix);
  lib_collection.Read(&read_pkg);
  timer.stop();
  xinfo("Read lib time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  auto insert_sizes = EstimateInsertSize(mapper, lib_collection);
  timer.stop();
  xinfo("Insert size estimation time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();

  MappingResultCollector collector(mapper.refseq().seq_count());
  MapToContigs(mapper, lib_collection, insert_sizes, &collector);
  timer.stop();
  xinfo("Mapping time elapsed: {}\n", timer.elapsed());

  timer.reset();
  timer.start();
  int32_t max_local_range = GetMaxLocalRange(lib_collection, insert_sizes);
  AssembleAndOutput(mapper, read_pkg, collector, opt.output_file,
                    max_local_range, opt);
  timer.stop();
  xinfo("Local assembly time elapsed: {}\n", timer.elapsed());
}
