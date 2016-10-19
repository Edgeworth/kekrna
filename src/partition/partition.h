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

#include <cmath>
#include "common.h"
#include "array.h"

namespace kekrna {
namespace partition {

enum : int8_t {
  PT_P,
  PT_U,
  PT_U2,
  PT_U_WC,
  PT_U_GU,
  PT_U_RCOAX,
  PT_SIZE
};

enum : int8_t {
  PTEXT,
  PTEXT_WC,     // Must start with a branch not involved in an interaction that is Watson-Crick
  PTEXT_GU,     // Must start with a branch not involved in an interaction that is GU
  PTEXT_RCOAX,  // Must start with a branch, that branch is involved backwards in a RCOAX stack.
  PTEXT_SIZE
};

// TODO rename this probability or something
typedef array3d_t<penergy_t, 1> probabilities_t;

struct partition_t {
  array3d_t<penergy_t, 1> p;
  penergy_t q;

//  partition_t(partition_t&&) = default; TODO
//  partition_t& operator=(partition_t&&) = default;
//
//  partition_t(const partition_t&) = delete;
//  partition_t& operator=(const partition_t&) = delete;
};

probabilities_t ComputeProbabilities(const partition_t& partition);

namespace internal {

void Partition0();
void Exterior();

}
}
}

#endif  // KEKRNA_PARTITION_H
