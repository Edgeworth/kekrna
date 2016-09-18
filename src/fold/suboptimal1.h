#ifndef KEKRNA_SUBOPTIMAL1_H
#define KEKRNA_SUBOPTIMAL1_H

#include <set>
#include <algorithm>
#include "common.h"
#include "parsing.h"
#include "fold/fold_internal.h"
#include "fold/suboptimal1_delta.h"
#include "fold/suboptimal1_num.h"
#include "energy/structure.h"

namespace kekrna {
namespace fold {
namespace internal {

class Suboptimal1 {
public:
  Suboptimal1(energy_t delta_, int num_) : delta(delta_), num(num_) {}

  std::vector<computed_t> Run() {
    verify_expr(gr.size() < std::numeric_limits<int16_t>::max(), "RNA too long for suboptimal folding");
    if (num == -1)
      return Suboptimal1Delta(delta).Run();
    return Suboptimal1Num(delta, num).Run();
  }
private:
  const energy_t delta;
  const int num;
};

}
}
}

#endif   // KEKRNA_SUBOPTIMAL1_H
