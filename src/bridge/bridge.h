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
  virtual std::vector<computed_t> Suboptimal(const primary_t& r, energy_t energy_delta) const = 0;
};

const std::map<std::string, ArgParse::option_t> BRIDGE_OPTIONS = {
    {"r", {"rnastructure"}}, {"k", {"kekrna"}}};

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse);
}
}

#endif  // KEKRNA_BRIDGE_H
