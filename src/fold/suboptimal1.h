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
#ifndef KEKRNA_SUBOPTIMAL1_H
#define KEKRNA_SUBOPTIMAL1_H

#include <algorithm>
#include "common.h"
#include "fold/fold.h"
#include "splaymap.h"

namespace kekrna {
namespace fold {
namespace internal {

struct expand_t {
  expand_t() = delete;
  expand_t(energy_t energy_) : energy(energy_) {}
  expand_t(energy_t energy_, const index_t& to_expand_) : energy(energy_), to_expand(to_expand_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const index_t& unexpanded_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const ctd_idx_t& ctd0_)
      : energy(energy_), to_expand(to_expand_), ctd0(ctd0_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const index_t& unexpanded_,
      const ctd_idx_t& ctd0_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_) {}
  expand_t(energy_t energy_, const index_t& to_expand_, const index_t& unexpanded_,
      const ctd_idx_t& ctd0_, const ctd_idx_t& ctd1_)
      : energy(energy_), to_expand(to_expand_), unexpanded(unexpanded_), ctd0(ctd0_), ctd1(ctd1_) {}
  energy_t energy;
  index_t to_expand, unexpanded;  // st is -1 if this does not exist
  ctd_idx_t ctd0, ctd1;

  bool operator<(const expand_t& o) const { return energy < o.energy; }
};

std::vector<expand_t> GenerateExpansions(const index_t& to_expand, energy_t delta);

class Suboptimal1 {
public:
  Suboptimal1(energy_t delta_, int num) :
      delta(delta_ == -1 ? CAP_E : delta_),
      max_structures(num == -1 ? MAX_STRUCTURES : num) {}

  int Run(SuboptimalCallback fn, bool sorted);

private:
  constexpr static int MAX_STRUCTURES = std::numeric_limits<int>::max() / 4;

  struct dfs_state_t {
    int idx;
    index_t expand;
    // Stores whether this node's |expand| was from |unexpanded| and needs to be replaced.
    bool should_unexpand;
  };

  const energy_t delta;
  const int max_structures;
  // This node is where we build intermediate results to be pushed onto the queue.
  SplayMap<index_t, std::vector<expand_t>> cache;
  std::vector<dfs_state_t> q;
  std::vector<index_t> unexpanded;

  std::pair<int, int> RunInternal(SuboptimalCallback fn,
      energy_t cur_delta, bool exact_energy, int structure_limit);

  const std::vector<expand_t>& GetExpansion(const index_t& to_expand) {
    if (!cache.Find(to_expand)) {
      // Need to generate the full way to delta so we can properly set |next_seen|.
      auto exps = GenerateExpansions(to_expand, delta);
      std::sort(exps.begin(), exps.end());
      auto res = cache.Insert(to_expand, std::move(exps));
      assert(res);
    }
    return cache.Get();
  }
};
}
}
}

#endif  // KEKRNA_SUBOPTIMAL1_H
