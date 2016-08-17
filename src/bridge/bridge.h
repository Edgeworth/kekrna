#ifndef KEKRNA_BRIDGE_H
#define KEKRNA_BRIDGE_H

#include "common.h"

#include "nn_unpaired_model.hpp"
#include "RNAstructure/rna_library.h"

namespace kekrna {
namespace bridge {

class RnaPackage {
public:

  virtual ~RnaPackage() = default;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const = 0;

  virtual folded_rna_t Fold(const rna_t& rna) const = 0;
};

class Rnark : public RnaPackage {
public:
  Rnark(const std::string& data_path);

  Rnark(const Rnark&) = delete;

  Rnark& operator=(const Rnark&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;

  virtual folded_rna_t Fold(const rna_t& rna) const;

private:
  librnary::NNUnpairedModel model;
};

class Rnastructure : public RnaPackage {
public:
  Rnastructure(const std::string& data_path, bool _use_lyngso);

  Rnastructure(const Rnastructure&) = delete;

  Rnastructure& operator=(Rnastructure&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;

  virtual folded_rna_t Fold(const rna_t& rna) const;

private:
  std::unique_ptr<datatable> data;
  bool use_lyngso;
};

class Kekrna : public RnaPackage {
public:
  Kekrna(const std::string& data_path);

  Kekrna(const Kekrna&) = delete;

  Kekrna& operator=(const Kekrna&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;

  virtual folded_rna_t Fold(const rna_t& rna) const;
};

}
}

#endif //KEKRNA_BRIDGE_H
