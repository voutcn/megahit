//
// Created by vout on 11/21/18.
//

#include "contig_writer.h"
#include "unitig_graph.h"
#include "definitions.h"
#include <cassert>

namespace {

static inline char Complement(char c) {
  if (c >= 0 && c < 4) {
    return 3 - c;
  }
  switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: assert(false);
  }
  return 0;
}

static inline void ReverseComplement(std::string &s) {
  int i, j;
  for (i = 0, j = s.length() - 1; i < j; ++i, --j) {
    std::swap(s[i], s[j]);
    s[i] = Complement(s[i]);
    s[j] = Complement(s[j]);
  }
  if (i == j) {
    s[i] = Complement(s[i]);
  }
}

void FoldPalindrome(std::string &s, unsigned kmer_k, bool is_loop) {
  if (is_loop) {
    for (unsigned i = 1; i + kmer_k <= s.length(); ++i) {
      std::string rc = s.substr(i, kmer_k);
      ReverseComplement(rc);
      if (rc == s.substr(i - 1, kmer_k)) {
        assert(i <= s.length() / 2);
        s = s.substr(i, s.length() / 2);
        break;
      }
    }
  } else {
    int num_edges = s.length() - kmer_k;
    assert(num_edges % 2 == 1);
    s.resize((num_edges - 1) / 2 + kmer_k + 1);
  }
}

}

void OutputContigs(UnitigGraph &graph, FILE *contig_file, FILE *final_file,
                   bool change_only, uint32_t min_standalone) {
  assert(!(change_only && final_file != nullptr)); // if output changed contigs, must not output final contigs

#pragma omp parallel for
  for (UnitigGraph::size_type i = 0; i < graph.size(); ++i) {
    auto adapter = graph.MakeVertexAdapter(i);
    double multi = change_only ? 1 : std::min(static_cast<double>(kMaxMul), adapter.avg_depth());
    std::string label = graph.VertexToDNAString(adapter);
    if (change_only && !adapter.is_changed()) {
      continue;
    }

    if (adapter.is_loop()) {
      int flag = contig_flag::kLoop | contig_flag::kStandalone;
      FILE *out_file = contig_file;

      if (adapter.is_palindrome()) {
        FoldPalindrome(label, graph.k(), adapter.is_loop());
        flag = contig_flag::kStandalone;
      }

      if (final_file != nullptr) {
        if (label.length() < min_standalone) {
          continue;
        } else {
          out_file = final_file;
        }
      }
      WriteContig(label, graph.k(), i, flag, multi, out_file);
    } else {
      FILE *out_file = contig_file;
      int flag = 0;

      if (adapter.forsure_standalone() || (graph.InDegree(adapter) == 0 && graph.OutDegree(adapter) == 0)) {
        if (adapter.is_palindrome()) {
          FoldPalindrome(label, graph.k(), adapter.is_loop());
        }
        flag = contig_flag::kStandalone;
        if (final_file != nullptr) {
          if (label.length() < min_standalone) {
            continue;
          } else {
            out_file = final_file;
          }
        }
      }
      WriteContig(label, graph.k(), i, flag, multi, out_file);
    }
  }
}