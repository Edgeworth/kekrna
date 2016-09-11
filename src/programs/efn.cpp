#include <cstdio>
#include "parsing.h"
#include "energy/structure.h"
#include "energy/load_model.h"

using namespace kekrna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions({
      {"v", {"verbose"}}
  });
  argparse.ParseOrExit(argc, argv);
  const auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 2, "requires primary sequence and dot bracket");

  const auto em = energy::LoadEnergyModelFromArgParse(argparse);
  std::unique_ptr<energy::Structure> structure;
  if (parsing::IsCtdString(pos.back())) {
    const auto computed = parsing::ParseCtdComputed(pos.front(), pos.back());
    printf("Energy: %d\n", energy::ComputeEnergyWithCtds(computed, *em, false, &structure).energy);
  } else {
    const auto secondary = parsing::ParseDotBracketSecondary(pos.front(), pos.back());
    printf("Energy: %d\n", energy::ComputeEnergy(secondary, *em, &structure).energy);
  }

  if (argparse.HasFlag("v")) {
    const auto descs = structure->Description();
    for (const auto& desc : descs) {
      printf("%s\n", desc.c_str());
    }
  }
}
