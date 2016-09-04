#ifndef KEKRNA_BRUTE_FOLD_H
#define KEKRNA_BRUTE_FOLD_H

#include "common.h"
#include "energy/energy_model.h"

namespace kekrna {
namespace fold {

computed_t FoldBruteForce(const primary_t& r, const energy::EnergyModel& em);

}
}
#endif  // KEKRNA_BRUTE_FOLD_H
