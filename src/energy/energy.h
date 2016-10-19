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
#ifndef KEKRNA_ENERGY_H
#define KEKRNA_ENERGY_H

#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <utility>
#include "argparse.h"
#include "common.h"
#include "energy/energy_model.h"

namespace kekrna {
namespace energy {

class Structure;

// If (st, en) is not paired, treated as an exterior loop.
energy_t ComputeSubstructureEnergy(computed_t& computed, bool compute_ctds,
    int st, int en, const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
computed_t ComputeEnergy(
    const secondary_t& secondary, const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
computed_t ComputeEnergyWithCtds(const computed_t& computed, const EnergyModel& em,
    bool compute_ctds = false, std::unique_ptr<Structure>* s = nullptr);
}
}

#endif  // KEKRNA_ENERGY_H
