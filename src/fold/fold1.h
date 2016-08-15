#ifndef KEKRNA_FOLD1_H
#define KEKRNA_FOLD1_H

#include "common.h"
#include "fold.h"
#include "array.h"

namespace kekrna {
namespace fold {

// Note that this doesn't completely follow the loaded energy model.
array3d_t<energy_t, DP_SIZE> ComputeTables1();

}
}

#endif //KEKRNA_FOLD1_H
