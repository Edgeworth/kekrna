#ifndef KEKRNA_BRIDGE_H
#define KEKRNA_BRIDGE_H

#include "common.h"

#include "nn_unpaired_model.hpp"
#include "RNAstructure/rna_library.h"

namespace kekrna {
namespace bridge {


class RnaPackage {
public:
  struct results_t {
    energy_t energy;
    folded_rna_t frna;
    std::string desc;
  };

  virtual ~RnaPackage() = default;
  virtual results_t Efn(const folded_rna_t& frna, bool verbose) = 0;
  virtual results_t Fold(const rna_t& rna, bool verbose) = 0;
};

class Rnark : public RnaPackage {
public:
  Rnark(const std::string& data_path);
  Rnark(const Rnark&) = delete;
  Rnark& operator=(const Rnark&) = delete;

  results_t Efn(const folded_rna_t& frna, bool verbose);
  results_t Fold(const rna_t& rna, bool verbose);

private:
  librnary::NNUnpairedModel model;
};

class Rnastructure : public RnaPackage {
public:
  Rnastructure(const std::string& data_path);
  Rnastructure(const Rnastructure&) = delete;
  Rnastructure& operator=(results_t&) = delete;

  results_t Efn(const folded_rna_t& frna, bool verbose);
  results_t Fold(const rna_t& rna, bool verbose);
private:
  std::unique_ptr<datatable> data;
};

class Kekrna : public RnaPackage {
public:
  Kekrna(const std::string& data_path);
  Kekrna(const Kekrna&) = delete;
  Kekrna& operator=(const Kekrna&) = delete;

  results_t Efn(const folded_rna_t& frna, bool verbose);
  results_t Fold(const rna_t& rna, bool verbose);
};

}
}

#endif //KEKRNA_BRIDGE_H
