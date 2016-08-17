#ifndef KEKRNA_FOLD_H
#define KEKRNA_FOLD_H

#include "constants.h"
#include "common.h"
#include "energy/energy.h"

namespace kekrna {
namespace fold {

enum {
  DP_P,  // For the paired array.
  DP_U,  // For the unpaired array.
  DP_U2, // Contains at least two branches.
  DP_U_WC,  // Unpaired but must start with a branch not involved in a CTD interaction that is not GU.
  DP_U_GU,  // Unpaired but must start with a branch not involved in a CTD interaction that is GU.
  DP_U_RCOAX,  // Unpaired but must start with a branch involved in a right coaxial stack - includes energy for it.
  DP_SIZE
};

enum {
  EXT,
  EXT_WC,
  EXT_GU,
  EXT_RCOAX,
  EXT_SIZE
};

const energy_t MULTILOOP_A = 93;

const energy_t AUGUBRANCH[4][4] = {
    {-6, -6, -6, 5 - 6},
    {-6, -6, -6, -6},
    {-6, -6, -6, 5 - 6},
    {5 - 6, -6, 5 - 6, -6}
};

energy_t FastTwoLoop(int ost, int oen, int ist, int ien);

inline bool IsNotLonely(int st, int en) {
  return (en - st - 3 >= constants::HAIRPIN_MIN_SZ && CanPair(r[st + 1], r[en - 1])) ||
      (st > 0 && en < int(r.size() - 1) && CanPair(r[st - 1], r[en + 1]));
}

energy_t Fold();

inline energy_t Fold(const rna_t& rna) {
  SetRna(rna);
  return Fold();
}

energy_t Fold2();
energy_t Fold1();
energy_t FoldSlow();
energy_t FoldBruteForce();

}
}

#endif //KEKRNA_FOLD_H