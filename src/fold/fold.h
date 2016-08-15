#ifndef KEKRNA_FOLD_H
#define KEKRNA_FOLD_H

#include "common.h"
#include "energy/energy.h"

namespace kekrna {
namespace fold {

energy_t Fold();
inline energy_t Fold(const rna_t& rna) {
  SetRna(rna);
  return Fold();
}
energy_t FoldBruteForce();

}
}

#endif //KEKRNA_FOLD_H
