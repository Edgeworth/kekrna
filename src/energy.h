#ifndef KEKRNA_ENERGY_H
#define KEKRNA_ENERGY_H

#include <utility>
#include "common.h"
#include "base.h"

namespace kekrna {
namespace energy {

typedef int32_t energy_t;
const energy_t AUGU_PENALTY = 5;

// Destroys folded_rna_t.
energy_t ComputeEnergy(folded_rna_t& frna);

}
}

#endif //KEKRNA_ENERGY_H
