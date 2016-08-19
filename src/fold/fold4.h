#ifndef KEKRNA_FOLD4_H
#define KEKRNA_FOLD4_H

#include "common.h"
#include "fold/fold.h"
#include "array.h"

namespace kekrna {
namespace fold {

// Note that this doesn't completely follow the loaded energy model.
array3d_t<energy_t, DP_SIZE> ComputeTables4();

}
}

#endif //KEKRNA_FOLD4_H
