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
#ifndef KEKRNA_PARTITION_GLOBALS_H
#define KEKRNA_PARTITION_GLOBALS_H

#include "common.h"
#include "energy/energy_model.h"
#include "partition/partition.h"

namespace kekrna {
namespace partition {
namespace internal {

extern array3d_t<penergy_t, PT_SIZE> gpt;
extern array2d_t<penergy_t, PTEXT_SIZE> gptext;

}

void SetPartitionGlobalState(const primary_t& r, const energy::EnergyModel& em);

}
}

#endif  // KEKRNA_PARTITION_GLOBALS_H