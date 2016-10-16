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
#ifndef KEKRNA_PARTITION_H
#define KEKRNA_PARTITION_H

#include "common.h"
#include "array.h"

namespace kekrna {
namespace partition {

enum {
  PT_P,
  PT_SIZE
};

typedef double penergy_t;
// TODO rename this probability or something
typedef array3d_t<penergy_t, 1> partition_t;

namespace internal {

void Partition0();

partition_t ComputeProbabilities();

}


}
}

#endif  // KEKRNA_PARTITION_H
