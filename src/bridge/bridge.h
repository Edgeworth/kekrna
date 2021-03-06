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
#ifndef KEKRNA_BRIDGE_H
#define KEKRNA_BRIDGE_H

#include "argparse.h"
#include "common.h"
#include "fold/context.h"

namespace kekrna {
namespace bridge {

class RnaPackage {
public:
  virtual ~RnaPackage() = default;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const = 0;
  virtual computed_t Fold(const primary_t& r) const = 0;
  virtual int Suboptimal(fold::SuboptimalCallback fn,
      const primary_t& r, energy_t energy_delta) const = 0;
  virtual std::vector<computed_t> SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const = 0;
};

const std::map<std::string, ArgParse::option_t> BRIDGE_OPTIONS = {
    {"r", {"rnastructure"}}, {"k", {"kekrna"}}};

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse);
}
}

#endif  // KEKRNA_BRIDGE_H
