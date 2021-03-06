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
#include <cstdio>
#include "energy/load_model.h"
#include "fold/context.h"
#include "parsing.h"

using namespace kekrna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  const auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  fold::Context ctx(parsing::StringToPrimary(pos.front()),
      energy::LoadEnergyModelFromArgParse(argparse), fold::ContextOptionsFromArgParse(argparse));
  const auto computed = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n", computed.energy, parsing::PairsToDotBracket(computed.s.p).c_str(),
      parsing::ComputedToCtdString(computed).c_str());
}
