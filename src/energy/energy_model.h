#ifndef KEKRNA_ENERGY_MODEL_H
#define KEKRNA_ENERGY_MODEL_H

#include <cstdarg>
#include <cstring>
#include "common.h"
#include "argparse.h"

namespace kekrna {
namespace energy {

const std::map<std::string, ArgParse::option_t> ENERGY_OPTIONS = {
    {"seed", ArgParse::option_t("seed for random energy model for kekrna").Arg()},
    {"data-path", ArgParse::option_t("data path for given energy model for kekrna").Arg("data/")}
};

void LoadEnergyModelFromDataDir(const std::string& data_dir);
void LoadRandomEnergyModel(int seed);
void LoadEnergyModelFromArgParse(const ArgParse& argparse);
uint32_t EnergyModelChecksum();


}
}

#endif //KEKRNA_ENERGY_MODEL_H
