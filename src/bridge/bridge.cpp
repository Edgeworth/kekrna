#include <constants.h>
#include "energy/energy.h"
#include "fold/fold.h"
#include "parsing.h"
#include "bridge.h"
#include "energy/structure.h"

#include "nn_unpaired_folder.hpp"
#include "nn_scorer.hpp"

namespace kekrna {
namespace bridge {

Rnark::Rnark(const std::string& data_path) : model(data_path) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
  model.SetMLParams(93, -6, 0, 0, 999999);
}

energy_t Rnark::Efn(const folded_rna_t& frna, std::string* desc) const {
  auto rna = librnary::StringToPrimary(parsing::RnaToString(frna.r));
  auto ss_tree = librnary::SSTree(librnary::DotBracketToMatching(parsing::PairsToDotBracket(frna.p)));

  librnary::NNScorer<librnary::NNUnpairedModel> scorer(model);
  scorer.SetRNA(rna);
  energy_t energy;
  if (desc) {
    auto trace = scorer.TraceExterior(ss_tree.RootSurface());
    *desc = trace.Describe(' ', 0) + "\n";
    energy = trace.recursive_score;
  } else {
    energy = scorer.ScoreExterior(ss_tree.RootSurface());
  }
  return energy;
}

folded_rna_t Rnark::Fold(const rna_t& rna) const {
  librnary::NNUnpairedFolder folder(model);
  folder.SetMaxTwoLoop(constants::TWOLOOP_MAX_SZ);
  folder.SetLonelyPairs(true);
  folder.SetStacking(true);
  auto primary = librnary::StringToPrimary(parsing::RnaToString(rna));
  energy_t energy = folder.Fold(primary);
  auto pairs = parsing::DotBracketToPairs(librnary::MatchingToDotBracket(folder.Traceback()));
  return {rna, pairs, energy};
}

Rnastructure::Rnastructure(const std::string& data_path, bool use_lyngso_) :
    data(librnary::LoadDatatable(data_path)), use_lyngso(use_lyngso_) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
}

energy_t Rnastructure::Efn(const folded_rna_t& frna, std::string* desc) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::RnaToString(frna.r)),
      librnary::DotBracketToMatching(parsing::PairsToDotBracket(frna.p))
  );
  efn2(data.get(), structure.get(), 1, true);
  return energy_t(structure->GetEnergy(1));
}

folded_rna_t Rnastructure::Fold(const rna_t& rna) const {
  auto structure = librnary::LoadStructure(
      librnary::StringToPrimary(parsing::RnaToString(rna)));
  // First false here says also generate the folding itself (not just the MFE).
  // Second last parameter is whether to generate the mfe structure only -- i.e. just one.
  // Last parameter is whether to use Lyngso or not.
  // Add two to TWOLOOP_MAX_SZ because rnastructure bug.
  dynamic(structure.get(), data.get(), 1, 0, 0, nullptr, false, nullptr, constants::TWOLOOP_MAX_SZ + 2, true,
      !use_lyngso);
  auto pairs = parsing::DotBracketToPairs(
      librnary::MatchingToDotBracket(librnary::StructureToMatching(*structure)));
  return {rna, pairs, energy_t(structure->GetEnergy(1))};
}

std::unique_ptr<Kekrna> KekrnaFromArgParse(const ArgParse& argparse) {
  energy_t (*fold_alg)() = &fold::Fold;
  auto opt = argparse.GetOption("alg");
  if (opt == "slow")
    fold_alg = &fold::FoldSlow;
  else if (opt == "1")
    fold_alg = &fold::Fold1;
  else if (opt == "2")
    fold_alg = &fold::Fold2;
  else if (opt == "3")
    fold_alg = &fold::Fold3;
  else if (opt == "4")
    fold_alg = &fold::Fold4;
  else if (opt == "brute")
    fold_alg = &fold::FoldBruteForce;
  else
    verify_expr(false, "unknown fold option");
  return std::make_unique<Kekrna>(argparse.GetOption("data-path"), fold_alg);
}

Kekrna::Kekrna(const std::string& data_path, energy_t (*fold_alg_)()) : fold_alg(fold_alg_) {
  verify_expr(data_path.size() && data_path.back() == '/', "invalid data path");
  LoadEnergyModelFromDataDir(data_path);
}

energy_t Kekrna::Efn(const folded_rna_t& frna, std::string* desc) const {
  energy_t energy;
  if (desc) {
    std::unique_ptr<structure::Structure> structure;
    energy = energy::ComputeEnergy(frna, &structure);
    for (const auto& s : structure->Description()) {
      *desc += s;
      *desc += "\n";
    }
  } else {
    energy = energy::ComputeEnergy(frna);
  }

  return energy;
}

folded_rna_t Kekrna::Fold(const rna_t& rna) const {
  SetRna(rna);
  auto energy = fold_alg();
  return {r, p, energy};
}

}
}
