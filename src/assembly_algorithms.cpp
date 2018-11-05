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

#include "assembly_algorithms.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <parallel/algorithm>
#include <unordered_set>
#include <queue>

#include "kmlib/bitvector.h"
#include "utils.h"
#include "histgram.h"

using std::vector;
using std::string;
using std::unordered_set;
using std::queue;

namespace assembly_algorithms {

double SetMinDepth(SuccinctDBG &dbg) {
  Histgram<mul_t> hist;

#pragma omp parallel for
  for (uint64_t i = 0; i < dbg.size(); ++i) {
    if (dbg.IsValidEdge(i)) {
      hist.insert(dbg.EdgeMultiplicity(i));
    }
  }

  double cov = hist.FirstLocalMinimum();
  for (int repeat = 1; repeat <= 100; ++repeat) {
    hist.TrimLow((mul_t) roundf(cov));
    unsigned median = hist.median();
    double cov1 = sqrt(median);
    if (abs(cov - cov1) < 1e-2) {
      return cov;
    }
    cov = cov1;
  }

  xwarning("Cannot detect min depth: unconverged");
  return 1;
}

int64_t Trim(SuccinctDBG &dbg, int len, AtomicBitVector &ignored) {
  int64_t number_tips = 0;
  AtomicBitVector to_remove(dbg.size());

#pragma omp parallel for reduction(+:number_tips)
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (!ignored.get(id) && dbg.EdgeOutdegreeZero(id)) {
      vector<uint64_t> path = {id};
      int64_t prev = -1;
      int64_t cur = id;
      bool is_tip = false;

      for (int i = 1; i < len; ++i) {
        prev = dbg.UniquePrevEdge(cur);
        if (prev == -1) {
          is_tip = dbg.EdgeIndegreeZero(cur);
          break;
        } else if (dbg.UniqueNextEdge(prev) == -1) {
          is_tip = true;
          break;
        } else {
          path.push_back(prev);
          cur = prev;
        }
      }
      if (is_tip) {
        for (unsigned long i : path) {
          to_remove.set(i);
        }
        ++number_tips;
        ignored.set(id);
        ignored.set(path.back());
        if (prev != -1) {
          ignored.unset(prev);
        }
      }
    }
  }

#pragma omp parallel for reduction(+:number_tips)
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (!ignored.get(id) && dbg.EdgeIndegreeZero(id)) {
      vector<uint64_t> path = {id};
      int64_t next = -1;
      int64_t cur = id;
      bool is_tip = false;

      for (int i = 1; i < len; ++i) {
        next = dbg.UniqueNextEdge(cur);
        if (next == -1) {
          is_tip = dbg.EdgeOutdegreeZero(cur);
          break;
        } else if (dbg.UniquePrevEdge(next) == -1) {
          is_tip = true;
          break;
        } else {
          path.push_back(next);
          cur = next;
        }
      }
      if (is_tip) {
        for (unsigned long i : path) {
          to_remove.set(i);
        }
        ++number_tips;
        ignored.set(id);
        ignored.set(path.back());
        if (next != -1) {
          ignored.unset(next);
        }
      }
    }
  }

#pragma omp parallel for
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (to_remove.get(id)) {
      dbg.SetInvalidEdge(id);
    }
  }
  return number_tips;
}

int64_t RemoveTips(SuccinctDBG &dbg, int max_tip_len, int min_final_standalone) {
  int64_t number_tips = 0;
  SimpleTimer timer;
  AtomicBitVector ignored(dbg.size());

#pragma omp parallel for
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (!dbg.EdgeIndegreeZero(id) && !dbg.EdgeOutdegreeZero(id)) {
      ignored.set(id);
    }
  }

  for (int len = 2; len < max_tip_len; len *= 2) {
    xlog("Removing tips with length less than %d; ", len);
    timer.reset();
    timer.start();
    number_tips += Trim(dbg, len, ignored);
    timer.stop();
    xlog_ext("Accumulated tips removed: %lld; time elapsed: %.4f\n", (long long) number_tips, timer.elapsed());
  }

  xlog("Removing tips with length less than %d; ", max_tip_len);
  timer.reset();
  timer.start();
  number_tips += Trim(dbg, max_tip_len, ignored);
  timer.stop();
  xlog_ext("Accumulated tips removed: %lld; time elapsed: %.4f\n", (long long) number_tips, timer.elapsed());

  return number_tips;
}

} // namespace assembly_algorithms