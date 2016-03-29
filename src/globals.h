#ifndef KEKRNA_GLOBALS_H
#define KEKRNA_GLOBALS_H

#include <vector>
#include "base.h"
#include "energy.h"

namespace kekrna {

// Globals.
// A X Y A
extern energy::energy_t stacking_e[4][4][4][4];
extern rna_t r;
extern std::vector<int> p;

}

#endif //KEKRNA_GLOBALS_H
