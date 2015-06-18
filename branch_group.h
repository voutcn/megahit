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

#ifndef BRANCH_GROUP_H_
#define BRANCH_GROUP_H_

#include <stdint.h>
#include <vector>
#include "succinct_dbg.h"

class BranchGroup {
  public:
    typedef std::vector<int64_t> BranchRecord;

    BranchGroup(SuccinctDBG *sdbg, int64_t begin_node, int max_branches = 2, int max_length = 0):
        sdbg_(sdbg), begin_node_(begin_node), max_branches_(max_branches), max_length_(max_length), status_(kNotMergedOrCorrected) {
        branches_.reserve(max_branches_);
        multiplicities_.reserve(max_branches_);
        if (max_length <= 0) {
            max_length_ = sdbg->kmer_k * 2 + 2;
        }
    }

    BranchGroup(const BranchGroup &rhs):
        sdbg_(rhs.sdbg_), begin_node_(rhs.begin_node_), end_node_(rhs.end_node_),
        max_branches_(rhs.max_branches_), max_length_(rhs.max_length_), branches_(rhs.branches_),
        multiplicities_(rhs.multiplicities_), status_(kNotMergedOrCorrected) { }

    BranchGroup& operator= (const BranchGroup &rhs) {
        sdbg_ = rhs.sdbg_;
        begin_node_ = rhs.begin_node_;
        end_node_ = rhs.end_node_;
        max_branches_ = rhs.max_branches_;
        max_length_ = rhs.max_length_;
        branches_ = rhs.branches_;
        multiplicities_ = rhs.multiplicities_;
        status_ = rhs.status_;

        return *this;
    }

    bool Search();
    bool RevSearch();
    bool RemoveErrorBranches(double cutoff_ratio = 0.5);
    bool Merge();
    size_t length() {
        if (branches_.size() == 0) {
            return 0;
        }
        return branches_[0].size();
    }

  private:
    SuccinctDBG *sdbg_;
    int64_t begin_node_;
    int64_t end_node_;
    int64_t max_branches_;
    int64_t max_length_;
    vector<BranchRecord> branches_;
    vector<int> multiplicities_;
    enum BranchGroupStatus {
        kNotMergedOrCorrected,
        kErrorRemoved,
        kCannotIdentifyError,
    } status_;
};

#endif