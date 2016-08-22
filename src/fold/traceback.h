#ifndef KEKRNA_TRACEBACK_H
#define KEKRNA_TRACEBACK_H

#include <stack>
#include "common.h"
#include "array.h"
#include "fold/fold.h"

namespace kekrna {
namespace fold {

typedef std::stack<std::tuple<int, int, int>> traceback_stack_t;

energy_t TraceExterior(const array3d_t<energy_t, DP_SIZE>& arr, traceback_stack_t& q);
void TraceStructure(const array3d_t<energy_t, DP_SIZE>& arr, traceback_stack_t& q);

}
}

#endif //KEKRNA_TRACEBACK_H
