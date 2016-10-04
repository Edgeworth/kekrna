// Copyright 2016, Eliot Courtney.
//
// This file is part of kekrna.
//
// kekrna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// kekrna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with kekrna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef KEKRNA_KEKRNA_H
#define KEKRNA_KEKRNA_H

#include "bridge/bridge.h"
#include "common.h"

namespace kekrna {
namespace bridge {

// Note that only one energy model can be loaded at a time.
class Kekrna : public RnaPackage {
public:
  Kekrna(energy::EnergyModelPtr em_, fold::context_options_t options_)
      : em(em_), options(options_) {}

  Kekrna(const Kekrna&) = delete;
  Kekrna& operator=(const Kekrna&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const override;
  virtual computed_t Fold(const primary_t& r) const override;
  virtual int Suboptimal(fold::SuboptimalCallback fn,
      const primary_t& r, energy_t energy_delta) const override;
  virtual std::vector<computed_t> SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const override;

private:
  const energy::EnergyModelPtr em;
  const fold::context_options_t options;
};
}
}

#endif  // KEKRNA_KEKRNA_H
