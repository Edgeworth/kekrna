#ifndef KEKRNA_SLOW_FOLD_H
#define KEKRNA_SLOW_FOLD_H

#include "common.h"
#include "array.h"

namespace kekrna {
namespace fold {

array3d_t<energy_t, DP_SIZE> ComputeTables_Slow();

}
}

#endif //KEKRNA_SLOW_FOLD_H
