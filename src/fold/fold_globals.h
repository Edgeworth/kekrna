#ifndef KEKRNA_FOLD_GLOBALS_H
#define KEKRNA_FOLD_GLOBALS_H

#include "common.h"

namespace kekrna {
namespace fold {

// For checks to make sure everything has been initialised.
extern bool g_fold_init;
extern energy_t g_augubranch[4][4];

}
}

#endif //KEKRNA_FOLD_GLOBALS_H
