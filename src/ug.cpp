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

#include "ug.h"

#include <omp.h>
#include <algorithm>
#include <string>
#include <utility>
#include <assembly/unitig_graph.h>
#include <assembly/contig_writer.h>
#include <assembly/bubble_remover.h>

#include "utils.h"
#include "definitions.h"
#include "sdbg/sdbg.h"
#include "kmlib/kmbitvector.h"

// -- helper functions --
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

// -- end of helper functions --

uint32_t UG::RemoveTips(uint32_t max_tip_len) {
  uint32_t num_removed = 0;
  for (uint32_t thre = 2; thre < max_tip_len; thre = std::min(thre * 2, max_tip_len)) {
#pragma omp parallel for schedule(static) reduction(+: num_removed)
    for (size_type i = 0; i < size(); ++i) {
      auto adapter = MakeVertexAdapter(i);
      if (adapter.length() >= thre) {
        continue;
      }
      VertexAdapter nexts[4], prevs[4];
      int outd = GetNextAdapters(adapter, nexts);
      int ind = GetPrevAdapters(adapter, prevs);

      if (ind + outd == 0) {
        adapter.MarkToDelete();
      } else if (outd == 1 && ind == 0) {
        if (nexts[0].avg_depth() > 8 * adapter.avg_depth()) {
          adapter.MarkToDelete();
        }
      } else if (outd == 0 && ind == 1) {
        if (prevs[0].avg_depth() > 8 * adapter.avg_depth()) {
          adapter.MarkToDelete();
        }
      }
      num_removed += adapter.is_to_delete();
    }

    Refresh(false);
    if (thre >= max_tip_len) { break; }
  }
  return num_removed;
}

int64_t UG::RemoveLowDepth(double min_depth) {
  int64_t num_removed = 0;
#pragma omp parallel for reduction(+: num_removed)
  for (size_type i = 0; i < size(); ++i) {
    auto adapter = MakeVertexAdapter(i);
    if (adapter.avg_depth() < min_depth) {
      adapter.MarkToDelete();
      ++num_removed;
    }
  }
  Refresh(false);
  return num_removed;
}

bool UG::RemoveLocalLowDepth(double min_depth,
                             uint32_t min_len,
                             uint32_t local_width,
                             double local_ratio,
                             int64_t &num_removed,
                             bool permanent_rm) {
  bool is_changed = false;
  bool need_refresh = false;

#pragma omp parallel for reduction(+: num_removed)
  for (size_type i = 0; i < size(); ++i) {
    auto adapter = MakeVertexAdapter(i);
    if (adapter.length() >= min_len) {
      continue;
    }
    int indegree = InDegree(adapter);
    int outdegree = OutDegree(adapter);
    if (indegree + outdegree == 0) {
      continue;
    }

    if ((indegree <= 1 && outdegree <= 1) || indegree == 0 || outdegree == 0) {
      double depth = adapter.avg_depth();
      if (is_changed && depth > min_depth)
        continue;
      double mean = LocalDepth(adapter, local_width);
      double threshold = min_depth;

      if (min_depth < mean * local_ratio)
        is_changed = true;
      else
        threshold = mean * local_ratio;

      if (depth < threshold) {
        is_changed = true;
        need_refresh = true;
        adapter.MarkToDelete();
        ++num_removed;
      }
    }
  }

  if (need_refresh) {
    bool set_changed = !permanent_rm;
    Refresh(set_changed);
  }
  return is_changed;
}

uint32_t UG::DisconnectWeakLinks(double local_ratio = 0.1) {
  uint32_t num_disconnected = 0;
#pragma omp parallel for reduction(+: num_disconnected)
  for (size_type i = 0; i < size(); ++i) {
    auto adapter = MakeVertexAdapter(i);
    if (adapter.is_loop() || adapter.is_palindrome()) {
      continue;
    }
    for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
      VertexAdapter nexts[4];
      double depths[4];
      double total_depth = 0;
      int degree = GetNextAdapters(adapter, nexts);
      if (degree <= 1) {
        continue;
      }
      for (int j = 0; j < degree; ++j) {
        depths[j] = nexts[j].avg_depth();
        total_depth += depths[j];
      }
      for (int j = 0; j < degree; ++j) {
        if (depths[j] <= local_ratio * total_depth) {
          num_disconnected += !nexts[j].is_to_disconnect();
          nexts[j].MarkToDisconnect();
        }
      }
    }
  }
  Refresh(false);
  return num_disconnected;
}

double UG::LocalDepth(VertexAdapter &adapter, uint32_t local_width) {
  double total_depth = 0;
  uint64_t num_added_edges = 0;

  for (int strand = 0; strand < 2; ++strand, adapter.ReverseComplement()) {
    VertexAdapter outs[4];
    int degree = GetNextAdapters(adapter, outs);

    for (int i = 0; i < degree; ++i) {
      if (outs[i].length() <= local_width) {
        num_added_edges += outs[i].length();
        total_depth += outs[i].total_depth();
      } else {
        num_added_edges += local_width;
        total_depth += outs[i].avg_depth() * local_width;
      }
    }
  }

  if (num_added_edges == 0) {
    return 0;
  } else {
    return total_depth / num_added_edges;
  }
}

void UG::OutputContigs(FILE *contig_file, FILE *final_file, Histgram<int64_t> &histo,
                       bool change_only, uint32_t min_final_standalone) {
  histo.clear();
  assert(!(change_only && final_file != nullptr)); // if output changed contigs, must not output final contigs

#pragma omp parallel for
  for (size_type i = 0; i < size(); ++i) {
    auto adapter = MakeVertexAdapter(i);
    double multi = std::min((double) kMaxMul, adapter.avg_depth());
    if (change_only) {
      multi = 1;
    }
    std::string label = VertexToDNAString(adapter);
    if (adapter.is_palindrome() && adapter.is_loop()) {
      FoldPalindrome(label, k(), adapter.is_loop());
    }
    histo.insert(label.length());

    if (change_only && !adapter.is_changed()) {
      continue;
    }

    if (adapter.is_loop()) {
//      adapter.MarkToDelete();
      int flag = contig_flag::kLoop | contig_flag::kIsolated;
      FILE *out_file = contig_file;

      if (adapter.is_palindrome()) {
        flag = contig_flag::kIsolated;
      }

      if (final_file != nullptr) {
        if (label.length() < min_final_standalone) {
          continue;
        } else {
          out_file = final_file;
        }
      }
      WriteContig(label, k(), i, flag, multi, out_file);
    } else {
      FILE *out_file = contig_file;
      int flag = 0;

      if (InDegree(adapter) == 0 && OutDegree(adapter) == 0) {
        adapter.MarkToDelete();
        if (adapter.is_palindrome()) {
          FoldPalindrome(label, k(), adapter.is_loop());
        }
        flag = contig_flag::kIsolated;
        if (final_file != nullptr) {
          if (label.length() < min_final_standalone) {
            continue;
          } else {
            out_file = final_file;
          }
        }
      }
      WriteContig(label, k(), i, flag, multi, out_file);
    }
  }
}

