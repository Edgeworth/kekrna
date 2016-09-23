#ifndef KEKRNA_KEKRNA_H
#define KEKRNA_KEKRNA_H

#include "common.h"
#include "bridge/bridge.h"

namespace kekrna {
namespace bridge {

// Note that only one energy model can be loaded at a time.
class Kekrna : public RnaPackage {
public:
  Kekrna(energy::EnergyModelPtr em_, fold::context_options_t options_) : em(em_), options(options_) {}

  Kekrna(const Kekrna&) = delete;
  Kekrna& operator=(const Kekrna&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const;
  virtual computed_t Fold(const primary_t& r) const;
  virtual std::vector<computed_t> Suboptimal(const primary_t& r, energy_t energy_delta) const;
private:
  const energy::EnergyModelPtr em;
  const fold::context_options_t options;
};


}
}

#endif  // KEKRNA_KEKRNA_H
