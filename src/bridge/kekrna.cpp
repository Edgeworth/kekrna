#include "bridge/kekrna.h"
#include "energy/structure.h"

namespace kekrna {
namespace bridge {

energy_t Kekrna::Efn(const secondary_t& secondary, std::string* desc) const {
  computed_t computed;
  if (desc) {
    std::unique_ptr<energy::Structure> structure;
    computed = energy::ComputeEnergy(secondary, *em, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    computed = energy::ComputeEnergy(secondary, *em);
  }

  return computed.energy;
}

computed_t Kekrna::Fold(const primary_t& r) const { return fold::Context(r, em, options).Fold(); }

std::vector<computed_t> Kekrna::Suboptimal(const primary_t& r, energy_t energy_delta) const {
  return fold::Context(r, em, options).Suboptimal(energy_delta, -1);
}
}
}
