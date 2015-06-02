/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
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

#include <assert.h>
#include "branch_group.h"

bool BranchGroup::Search() {
    if (sdbg_->Indegree(begin_node_) != 1 || sdbg_->Outdegree(begin_node_) > max_branches_ || sdbg_->Outdegree(begin_node_) <= 1) {
        return false;
    }

    branches_.push_back(BranchRecord(1, begin_node_));
    multiplicities_.push_back(0);

    bool converged = false;

    for (int j = 1; j < max_length_; ++j) {
        int num_branches = branches_.size();
        for (int i = 0; i < num_branches; ++i) {
            int64_t current = branches_[i].back();
            int64_t outgoings[4];
            int edge_countings[4];
            int out_degree = sdbg_->Outgoings(current, outgoings, edge_countings);

            for (int x = 0; x < out_degree; ++x) {
                assert(sdbg_->IsLast(outgoings[x]));
            }

            if (out_degree >= 1) {
                // append this node the last of current branch
                branches_[i].push_back(outgoings[0]);
                multiplicities_[i] += edge_countings[0];

                if ((int)branches_.size() + out_degree - 1 > max_branches_) {
                    // too many branches
                    return false;
                } else {
                    BranchRecord curr_branch = branches_[i];
                    int curr_branch_multiplicity = multiplicities_[i] - edge_countings[0];
                    for (int x = 1; x < out_degree; ++x) {
                        curr_branch.pop_back();
                        curr_branch.push_back(outgoings[x]);
                        branches_.push_back(curr_branch);
                        multiplicities_.push_back(curr_branch_multiplicity + edge_countings[x]);
                    }
                }
            }
        }

        // check whether all branches's last nodes are coming from this branch group
        for (unsigned i = 0; i < branches_.size(); ++i) {
            int64_t last_node = branches_[i].back();
            int64_t incomings[4];
            int in_degree = sdbg_->Incomings(last_node, incomings);
            for (int x = 0; x < in_degree; ++x) {
                assert(sdbg_->IsLast(incomings[x]));
            }

            if (in_degree == 1) {
                continue;
            } else {
                for (int x = 0; x < in_degree; ++x) {
                    bool exist_in_group = false;
                    for (auto it = branches_.begin(); it != branches_.end(); ++it) {
                        if ((*it)[j - 1] == incomings[x]) {
                            exist_in_group = true;
                            break;
                        }
                    }
                    if (!exist_in_group) {
                        return false;
                    }
                }
            }
        }

        // check converge
        end_node_ = branches_[0].back();
        if (sdbg_->Outdegree(end_node_) == 1) {
            converged = true;
            for (unsigned i = 1; i < branches_.size(); ++i) {
                if (branches_[i].back() != end_node_) {
                    converged = false;
                    break;
                }
            }
            if (converged) {
                break;
            }
        }
    }
    return converged && begin_node_ != end_node_;
}

bool BranchGroup::RemoveErrorBranches(double cutoff_ratio) {
    int best_multiplicity = multiplicities_[0];

    for (unsigned i = 1; i < branches_.size(); ++i) {
        int curr_multiplicity = multiplicities_[i];
        if (curr_multiplicity >= best_multiplicity) {
            best_multiplicity = curr_multiplicity;
        }
    }

    vector<unsigned> not_removed;
    not_removed.reserve(branches_.size());
    for (unsigned i = 0; i < branches_.size(); ++i) {
        if (multiplicities_[i] == best_multiplicity || multiplicities_[i] > cutoff_ratio * best_multiplicity) {
            not_removed.push_back(i);
            continue;
        }
        for (unsigned j = 1; j + 1 < branches_[i].size(); ++j) {
            sdbg_->SetInvalid(branches_[i][j]);
        }
    }

/*
    if (not_removed.size() > 1) {
        int64_t biggest_id = -1;
        unsigned remain_branch = 0;
        for (unsigned i = 0; i < not_removed.size(); ++i) {
            int j = not_removed[i];
            int64_t b = -1;

            for (int x = 1; x < branches_[j].size() - 1; ++x) {
                int64_t node = branches_[j][x];
                b = std::max(std::max(b, node), sdbg_->ReverseComplement(node));
            }

            if (b > biggest_id) {
                biggest_id = b;
                remain_branch = i;
            }
        }

        for (unsigned i = 0; i < not_removed.size(); ++i) {
            if (i == remain_branch) {
                continue;
            }
            for (unsigned j = 1; j + 1 < branches_[not_removed[i]].size(); ++j) {
                sdbg_->SetInvalid(branches_[not_removed[i]][j]);
            }
        }

        not_removed[0] = not_removed[remain_branch];
        not_removed.resize(1);        
    }
*/
    unsigned num_remained = 0;
    for (unsigned i = 0; i < not_removed.size(); ++i) {
        if (num_remained != not_removed[i]) {
            branches_[num_remained] = branches_[not_removed[i]];
        }
        for (unsigned j = 1; j + 1 < branches_[num_remained].size(); ++j) {
            sdbg_->SetValid(branches_[num_remained][j]);
        }
        ++num_remained;
    }

    branches_.resize(num_remained);

    if (num_remained == 1) {
        status_ = kErrorRemoved;
        return true;
    } else {
        status_ = kCannotIdentifyError;
        return false;
    }
}