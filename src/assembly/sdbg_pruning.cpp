/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics
 * Limited
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

#include "sdbg_pruning.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <vector>

#include "kmlib/kmbitvector.h"
#include "utils/histgram.h"
#include "utils/utils.h"

namespace sdbg_pruning {

double InferMinDepth(SDBG &dbg) {
  Histgram<mul_t> hist;

#pragma omp parallel for
  for (uint64_t i = 0; i < dbg.size(); ++i) {
    if (dbg.IsValidEdge(i)) {
      hist.insert(dbg.EdgeMultiplicity(i));
    }
  }

  double cov = hist.FirstLocalMinimum();
  for (int repeat = 1; repeat <= 100; ++repeat) {
    hist.TrimLow(static_cast<mul_t>(roundf(cov)));
    unsigned median = hist.median();
    double cov1 = sqrt(median);
    if (fabs(cov - cov1) < 1e-2) {
      return cov;
    }
    cov = cov1;
  }

  xwarn("Cannot detect min depth: unconverged");
  return 1;
}

int64_t Trim(SDBG &dbg, int len, AtomicBitVector &ignored) {
  int64_t number_tips = 0;
  AtomicBitVector to_remove(dbg.size());
  std::vector<uint64_t> path;

#pragma omp parallel for reduction(+ : number_tips) private(path)
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (!ignored.at(id) && dbg.EdgeOutdegreeZero(id)) {
      uint64_t prev = SDBG::kNullID;
      uint64_t cur = id;
      bool is_tip = false;
      path.clear();
      path.push_back(id);

      for (int i = 1; i < len; ++i) {
        prev = dbg.UniquePrevEdge(cur);
        if (prev == SDBG::kNullID) {
          is_tip = dbg.EdgeIndegreeZero(cur);
          break;
        } else if (dbg.UniqueNextEdge(prev) == SDBG::kNullID) {
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
        if (prev != SDBG::kNullID) {
          ignored.unset(prev);
        }
      }
    }
  }

#pragma omp parallel for reduction(+ : number_tips) private(path)
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (!ignored.at(id) && dbg.EdgeIndegreeZero(id)) {
      uint64_t next = SDBG::kNullID;
      uint64_t cur = id;
      bool is_tip = false;
      path.clear();
      path.push_back(id);

      for (int i = 1; i < len; ++i) {
        next = dbg.UniqueNextEdge(cur);
        if (next == SDBG::kNullID) {
          is_tip = dbg.EdgeOutdegreeZero(cur);
          break;
        } else if (dbg.UniquePrevEdge(next) == SDBG::kNullID) {
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
        if (next != SDBG::kNullID) {
          ignored.unset(next);
        }
      }
    }
  }

#pragma omp parallel for
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (to_remove.at(id)) {
      dbg.SetInvalidEdge(id);
    }
  }
  return number_tips;
}

uint64_t RemoveTips(SDBG &dbg, int max_tip_len) {
  uint64_t number_tips = 0;
  SimpleTimer timer;
  AtomicBitVector ignored(dbg.size());

#pragma omp parallel for
  for (uint64_t id = 0; id < dbg.size(); ++id) {
    if (!dbg.EdgeIndegreeZero(id) && !dbg.EdgeOutdegreeZero(id)) {
      ignored.set(id);
    }
  }

  for (int len = 2; len < max_tip_len; len *= 2) {
    xinfo("Removing tips with length less than {}; ", len);
    timer.reset();
    timer.start();
    number_tips += Trim(dbg, len, ignored);
    timer.stop();
    xinfoc("Accumulated tips removed: {}; time elapsed: {.4}\n", number_tips,
           timer.elapsed());
  }

  xinfo("Removing tips with length less than {}; ", max_tip_len);
  timer.reset();
  timer.start();
  number_tips += Trim(dbg, max_tip_len, ignored);
  timer.stop();
  xinfoc("Accumulated tips removed: {}; time elapsed: {.4}\n", number_tips,
         timer.elapsed());

  return number_tips;
}

}  // namespace sdbg_pruning