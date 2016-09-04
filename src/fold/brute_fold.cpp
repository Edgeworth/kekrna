#include "fold/brute_fold.h"
#include "fold/fold.h"

namespace kekrna {
namespace fold {

namespace {

computed_t best_computed;
secondary_t cur_secondary;
std::vector<std::pair<int, int>> base_pairs;
const energy::EnergyModel* em;

void FoldBruteForceInternal(int idx) {
  if (idx == int(base_pairs.size())) {
    auto computed = energy::ComputeEnergy(cur_secondary, *em);
    if (computed.energy < best_computed.energy)
      best_computed = std::move(computed);
    return;
  }
  // Don't take this base pair.
  FoldBruteForceInternal(idx + 1);

  // Take this base pair.
  bool can_take = true;
  const auto& pair = base_pairs[idx];
  // Only need to check in the range of this base pair. Since we ordered by
  // increasing st, anything at or after this will either be the start of something starting at st, or
  // something ending, both of which conflict with this base pair.
  for (int i = pair.first; i <= pair.second; ++i) {
    if (cur_secondary.p[i] != -1) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    cur_secondary.p[pair.first] = pair.second;
    cur_secondary.p[pair.second] = pair.first;
    FoldBruteForceInternal(idx + 1);
    cur_secondary.p[pair.first] = -1;
    cur_secondary.p[pair.second] = -1;
  }
}

}

computed_t FoldBruteForce(const primary_t& r, const energy::EnergyModel& em_ref) {
  best_computed = computed_t(r);
  cur_secondary = secondary_t(r);
  base_pairs.clear();
  em = &em_ref;
  // Add base pairs in order of increasing st, then en.
  for (int st = 0; st < int(r.size()); ++st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < int(r.size()); ++en) {
      if (internal::ViableFoldingPair(r, st, en))
        base_pairs.emplace_back(st, en);
    }
  }
  FoldBruteForceInternal(0);
  return std::move(best_computed);
}

}
}
