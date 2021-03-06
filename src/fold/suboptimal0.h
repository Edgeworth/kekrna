// Copyright 2016, Eliot Courtney.
//
// This file is part of kekrna.
//
// kekrna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// kekrna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with kekrna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef KEKRNA_SUBOPTIMAL0_H
#define KEKRNA_SUBOPTIMAL0_H

#include <set>
#include "common.h"
#include "energy/structure.h"
#include "fold/fold.h"
#include "parsing.h"

namespace kekrna {
namespace fold {
namespace internal {

class Suboptimal0 {
public:
  Suboptimal0(energy_t delta_, int num)
      : max_energy(delta_ == -1 ? CAP_E : gext[0][EXT] + delta_),
        max_structures(num == -1 ? std::numeric_limits<int>::max() / 4 : num) {
    verify_expr(max_structures > 0, "must request at least one structure");
  }
  int Run(SuboptimalCallback fn);

private:
  struct node_t {
    // State should be fully defined by |not_yet_expanded|, |history|, and |base_ctds| which denote
    // what it has done so far, and what it can do from now.
    std::vector<index_t> not_yet_expanded;
    std::vector<index_t> history;
    std::vector<int16_t> p;
    std::vector<Ctd> base_ctds;
    energy_t energy;  // Stores the minimum energy this state could have.

    bool operator<(const node_t& o) const { return energy < o.energy; }
  };

  const energy_t max_energy;
  const int max_structures;
  // This node is where we build intermediate results to be pushed onto the queue.
  node_t curnode;
  std::multiset<node_t> finished;
  std::multiset<node_t> q;

  void PruneInsert(std::multiset<node_t>& prune, const node_t& node) {
    if (node.energy <= max_energy) {
      if (int(prune.size()) >= max_structures && (--prune.end())->energy > node.energy)
        prune.erase(--prune.end());
      if (int(prune.size()) < max_structures) prune.insert(node);
    }
  }

  // Creates and inserts a new node with energy |energy| that doesn't
  // need to expand any more ranges than it currently has.
  void Expand(energy_t energy) {
    curnode.energy = energy;
    PruneInsert(q, curnode);
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the given ranges.
  void Expand(energy_t energy, index_t nye) {
    curnode.not_yet_expanded.push_back(nye);
    curnode.energy = energy;
    PruneInsert(q, curnode);
    curnode.not_yet_expanded.pop_back();
  }

  void Expand(energy_t energy, index_t nye, ctd_idx_t ctd_idx) {
    curnode.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye);
    curnode.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  // Creates and inserts a new node with energy |energy| that needs to expand the two given ranges.
  void Expand(energy_t energy, index_t nye0, index_t nye1) {
    curnode.not_yet_expanded.push_back(nye0);
    curnode.not_yet_expanded.push_back(nye1);
    curnode.energy = energy;
    PruneInsert(q, curnode);
    curnode.not_yet_expanded.pop_back();
    curnode.not_yet_expanded.pop_back();
  }

  void Expand(energy_t energy, index_t nye0, index_t nye1, ctd_idx_t ctd_idx) {
    curnode.base_ctds[ctd_idx.idx] = ctd_idx.ctd;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx.idx] = CTD_NA;
  }

  void Expand(energy_t energy, index_t nye0, index_t nye1, ctd_idx_t ctd_idx0, ctd_idx_t ctd_idx1) {
    curnode.base_ctds[ctd_idx0.idx] = ctd_idx0.ctd;
    curnode.base_ctds[ctd_idx1.idx] = ctd_idx1.ctd;
    Expand(energy, nye0, nye1);
    curnode.base_ctds[ctd_idx0.idx] = CTD_NA;
    curnode.base_ctds[ctd_idx1.idx] = CTD_NA;
  }
};
}
}
}

#endif  // KEKRNA_SUBOPTIMAL0_H
