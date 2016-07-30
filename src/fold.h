#ifndef KEKRNA_FOLD_H
#define KEKRNA_FOLD_H

#include "common.h"
#include "energy.h"

namespace kekrna {
namespace fold {

// This function is not re-entrant.
energy_t Fold(const rna_t& rna, std::unique_ptr<structure::Structure>* s = nullptr);

}
}

#endif //KEKRNA_FOLD_H
