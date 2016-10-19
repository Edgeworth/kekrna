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
#include <algorithm>
#include "fold/brute_fold.h"
#include "bridge/rnastructure.h"
#include "energy/load_model.h"
#include "partition/partition_globals.h"
#include "parsing.h"

using namespace kekrna;

void PrintPartition(const partition::partition_t& p,
    const std::string& name) {
  const int N = int(p.p.Size());
  printf("%s total: %lf\n", name.c_str(), p.q);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      printf("%lf ", p.p[i][j][0]);
    }
    printf("\n");
  }
}

int main(int argc, char* argv[]) {
  ArgParse argparse({{"cutoff", opt_t("Only return probabilities larger than this.").Arg("0.1")}});
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.AddOptions(CONTEXT_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  const auto& pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");
  const auto primary = parsing::StringToPrimary(pos.front());
  const auto em = energy::LoadEnergyModelFromArgParse(argparse);
  const auto cutoff = atof(argparse.GetOption("cutoff").c_str());
  const bridge::Rnastructure rnastructure("extern/miles_rnastructure/data_tables/", false);

  PrintPartition(rnastructure.Partition(primary).first, "RNAstructure");

  Context ctx(parsing::StringToPrimary(pos.front()), em, ContextOptionsFromArgParse(argparse));
  PrintPartition(ctx.Partition(), "kekrna");

  PrintPartition(fold::PartitionBruteForce(primary, *em).first, "Brute force");

//  std::vector<std::tuple<double, int, int>> pairs;
//  for (int i = 0; i < int(primary.size()); ++i) {
//    for (int j = 0; j < int(primary.size()); ++j)
//      if (probs[i][j][0] >= cutoff)
//        pairs.emplace_back(probs[i][j][0], i, j);
//  }
//  std::sort(pairs.rbegin(), pairs.rend());
//  for (const auto& pair : pairs) {
//    int st, en;
//    double prob;
//    std::tie(prob, st, en) = pair;
//    printf("%d %d: %.4lf %.4lf\n", st, en, prob, -log10(prob));
//  }
}
