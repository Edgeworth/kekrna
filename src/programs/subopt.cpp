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
  argparse.AddOptions({
      {"delta", ArgParse::option_t("maximum energy delta from minimum").Arg("-1")},
      {"num", ArgParse::option_t("maximum number of reported structures").Arg("-1")},
      {"q", ArgParse::option_t("quiet")},
      {"sorted", ArgParse::option_t("if the structures should be sorted")},
      {"ctd-output", ArgParse::option_t("if we should output CTD data")}
  });
  argparse.ParseOrExit(argc, argv);
  const auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  auto opt = fold::ContextOptionsFromArgParse(argparse);
  opt.table_alg = fold::context_options_t::TableAlg::TWO;
  fold::Context ctx(
      parsing::StringToPrimary(pos.front()), energy::LoadEnergyModelFromArgParse(argparse), opt);

  const energy_t subopt_delta = atoi(argparse.GetOption("delta").c_str());
  const int subopt_num = atoi(argparse.GetOption("num").c_str());
  const bool should_print = !argparse.HasFlag("q");
  const bool sorted = argparse.HasFlag("sorted");
  const bool ctd_data = argparse.HasFlag("ctd-output");
  verify_expr(subopt_delta >= 0 || subopt_num > 0, "nothing to do");

  fold::SuboptimalCallback fn = [](const computed_t&) {};
  if (should_print) {
    if (ctd_data) {
      fn = [](const computed_t& c) {
        printf("%d ", c.energy);
        puts(parsing::ComputedToCtdString(c).c_str());
      };
    } else {
      fn = [](const computed_t& c) {
        printf("%d ", c.energy);
        puts(fold::internal::grep.c_str());  // kek
      };
    }
  }
  int num_structures = ctx.Suboptimal(fn, sorted, subopt_delta, subopt_num);
  printf("%d suboptimal structures\n", num_structures);
}
