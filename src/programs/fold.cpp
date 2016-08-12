#include <cstdio>
#include <cassert>
#include "base.h"
#include "parsing.h"
#include "structure.h"
#include "fold.h"

using namespace kekrna;

int main(int argc, char* argv[]) {
  verify_expr(argc == 2, "requires one argument");
  LoadEnergyModelFromDataDir();
  auto rna = parsing::ParseRnaFromString(argv[1]);
  energy_t e = fold::Fold(rna);
  printf("Energy: %d\n%s\n", e, parsing::DotBracketFromPairs(p).c_str());
}
