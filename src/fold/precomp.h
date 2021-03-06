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
#ifndef KEKRNA_PRECOMP_H
#define KEKRNA_PRECOMP_H

#include "common.h"
#include "energy/energy_model.h"

namespace kekrna {
namespace fold {
namespace internal {

struct hairpin_precomp_t {
  static const int MAX_SPECIAL_HAIRPIN_SZ = 6;
  hairpin_precomp_t() : num_c(0) { memset(special, MAX_E & 0xFF, sizeof(special)); }

  energy_t special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};

struct precomp_t {
  energy_t augubranch[4][4];
  energy_t min_mismatch_coax;
  energy_t min_flush_coax;
  energy_t min_twoloop_not_stack;

  std::vector<hairpin_precomp_t> hairpin;
};

int MaxNumContiguous(const primary_t& r);
precomp_t PrecomputeData(const primary_t& r, const energy::EnergyModel& em);
// Must have global state set.
energy_t FastTwoLoop(int ost, int oen, int ist, int ien);
energy_t FastHairpin(int st, int en);
}
}
}

#endif  // KEKRNA_PRECOMP_H
