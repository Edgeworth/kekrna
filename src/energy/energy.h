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

computed_t ComputeEnergy(
    const secondary_t& secondary, const EnergyModel& em, std::unique_ptr<Structure>* s = nullptr);
computed_t ComputeEnergyWithCtds(const computed_t& computed, const EnergyModel& em,
    bool compute_ctds = false, std::unique_ptr<Structure>* s = nullptr);
}
}

#endif  // KEKRNA_ENERGY_H
