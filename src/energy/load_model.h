#ifndef KEKRNA_LOAD_MODEL_H
#define KEKRNA_LOAD_MODEL_H

#include "common.h"
#include "energy/energy_model.h"

namespace kekrna {
namespace energy {

const std::map<std::string, ArgParse::option_t> ENERGY_OPTIONS = {
    {"seed", ArgParse::option_t("seed for random energy model for kekrna").Arg()},
    {"data-path", ArgParse::option_t("data path for given energy model for kekrna").Arg("data/")}
};

EnergyModel LoadEnergyModelFromDataDir(const std::string& data_dir);
EnergyModel LoadRandomEnergyModel(uint_fast32_t seed);
EnergyModel LoadEnergyModelFromArgParse(const ArgParse& argparse);

}
}

#endif   //KEKRNA_LOAD_MODEL_H
