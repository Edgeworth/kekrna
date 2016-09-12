#include <cstdio>
#include <random>
#include <chrono>
#include <cinttypes>
#include <set>
#include "fold/brute_fold.h"
#include "fold/globals.h"
#include "fold/context.h"
#include "bridge/bridge.h"
#include "parsing.h"
#include "energy/load_model.h"

using namespace kekrna;
using namespace fold;
using namespace fold::internal;

const int SUBOPT_BRUTE_MAX_STRUCTURES = 10000;
const energy_t SUBOPT_MAX_DELTA = 6;  // Same as RNAstructure default.

#include "energy/structure.h"

class Fuzzer {
public:
  typedef std::deque<std::string> error_t;

  Fuzzer(const primary_t& r_, const energy::EnergyModelPtr em_, bool random_model_, uint_fast32_t seed_,
      const std::vector<context_options_t>& kekrnas_, const bridge::Rnastructure& rnastructure_)
      : r(r_), em(random_model_ ? energy::LoadRandomEnergyModel(seed) : em_),
        random_model(random_model_), seed(seed_), kekrnas(kekrnas_), rnastructure(rnastructure_) {}

  error_t Run() {
    error_t errors;
    AppendErrors(errors, MaybePrependHeader(KekrnaComputeAndCheckState(), "kekrna:"));
    if (!random_model) {
      AppendErrors(errors, MaybePrependHeader(RnastructureComputeAndCheckState(), "rnastructure:"));

      // Suboptimal folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does strange things
      // when the energy for suboptimal structures is 0 or above.
      if (kekrna_computeds[0].energy < -SUBOPT_MAX_DELTA) {
        context_options_t options(context_options_t::TableAlg::TWO, context_options_t::SuboptimalAlg::ZERO);
        Context ctx(r, em, options);
        auto kekrna_subopt = ctx.Suboptimal(SUBOPT_MAX_DELTA, -1);
        const auto rnastructure_subopt = rnastructure.Suboptimal(r, SUBOPT_MAX_DELTA);
        AppendErrors(errors, MaybePrependHeader(CheckSuboptimal(kekrna_subopt, true), "kekrna suboptimal:"));
        AppendErrors(errors,
            MaybePrependHeader(CheckSuboptimal(rnastructure_subopt, false), "rnastructure suboptimal:"));
        AppendErrors(errors, MaybePrependHeader(CheckSuboptimalPair(
            kekrna_subopt, rnastructure_subopt), "kekrna vs rnastructure suboptimal:"));
      }
    }
    if (r.size() <= 25)
      AppendErrors(errors, MaybePrependHeader(CheckBruteForce(), "brute force:"));
    AppendErrors(errors, MaybePrependHeader(CheckDpTables(), "dp tables:"));

    if (!errors.empty()) {
      if (random_model)
        errors.push_front(sfmt("Used random energy model with seed: %" PRIuFAST32 "\n", seed));
      else
        errors.push_front(sfmt("Used T04 energy model"));
      errors = MaybePrependHeader(errors,
          sfmt("Difference on len %zu RNA %s:", r.size(), parsing::PrimaryToString(r).c_str()));
    }

    return errors;
  }

private:
  const primary_t& r;
  const energy::EnergyModelPtr em;
  const bool random_model;
  const uint_fast32_t seed;
  const std::vector<context_options_t>& kekrnas;
  const bridge::Rnastructure& rnastructure;

  // Fuzz state.
  std::vector<computed_t> kekrna_computeds;
  std::vector<array3d_t<energy_t, DP_SIZE>> kekrna_dps;
  dp_state_t rnastructure_dp;

  error_t MaybePrependHeader(const error_t& main, const std::string& header) {
    if (main.empty()) return main;
    error_t nmain;
    nmain.push_front(header);
    for (auto& error : main)
      nmain.push_back("  " + error);  // mfw this inefficiency
    return nmain;
  }

  void AppendErrors(error_t& main, error_t&& extra) {
    for (auto& s : extra)
      main.push_back(std::move(s));
  }

  bool HasDuplicates(const std::vector<computed_t>& computeds) {
    // If energies are different but everything else is the same, it is still a bug.
    std::set<std::pair<secondary_t, std::vector<Ctd>>> suboptimal_set;
    for (const auto& computed : computeds) {
      auto val = std::make_pair(computed.s, computed.base_ctds);
      if (suboptimal_set.count(val))
        return true;
      suboptimal_set.insert(val);
    }
    return false;
  }

  error_t CheckSuboptimal(const std::vector<computed_t>& subopt, bool has_ctds) {
    error_t errors;
    // Check at least one suboptimal structure.
    if (subopt.empty())
      errors.push_back("no structures returned");
    // Check MFE.
    if (!subopt.empty() && kekrna_computeds[0].energy != subopt[0].energy)
      errors.push_back(sfmt("lowest structure energy %d != mfe %d", subopt[0].energy, kekrna_computeds[0].energy));

    // Only ones with CTDs set can do these tests.
    if (has_ctds) {
      // Check for duplicate structures.
      if (HasDuplicates(subopt))
        errors.push_back("has duplicates");

      for (int i = 0; i < int(subopt.size()); ++i) {
        const auto& structure = subopt[i];
        auto suboptimal_efn = energy::ComputeEnergyWithCtds(structure, *em);
        if (suboptimal_efn.energy != structure.energy) {
          errors.push_back(sfmt("structure %d: energy %d != efn %d", i, structure.energy, suboptimal_efn.energy));
          break;
        }

        // Incidentally test ctd parsing.
        auto parsed_computed = parsing::ParseCtdComputed(
            parsing::PrimaryToString(structure.s.r), parsing::ComputedToCtdString(structure));
        parsed_computed.energy = structure.energy;
        if (parsed_computed != structure) {
          errors.push_back(sfmt("structure %d: bug in parsing code", i));
          break;
        }
      }
    }
    return errors;
  }

  error_t CheckSuboptimalPair(const std::vector<computed_t>& a, const std::vector<computed_t>& b) {
    error_t errors;
    if (a.size() != b.size()) {
      errors.push_back(sfmt("first has %d structures != second has %d structures", int(a.size()), int(b.size())));
    } else {
      for (int i = 0; i < int(a.size()); ++i) {
        if (a[i].energy != b[i].energy) {
          errors.push_back(sfmt("structure %d: first %d != second %d", i, a[i].energy, b[i].energy));
          break;
        }
      }
    }
    return errors;
  }

  error_t CheckDpTables() {
    error_t errors;
    int st = 0, en = 0, a = 0;
    const int N = int(r.size());
    for (st = N - 1; st >= 0; --st) {
      for (en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        for (a = 0; a < DP_SIZE; ++a) {
          const auto kekrna0 = kekrna_dps[0][st][en][a];
          for (int i = 0; i < int(kekrnas.size()); ++i) {
            const auto kekrnai = kekrna_dps[i][st][en][a];
            // If meant to be infinity and not.
            if (((kekrna0 < constants::CAP_E) != (kekrnai < constants::CAP_E)) ||
                (kekrna0 < constants::CAP_E && kekrna0 != kekrnai)) {
              errors.push_back(sfmt("kekrna %d at %d %d %d: %d != %d",
                  i, st, en, a, kekrna_dps[i][st][en][a], kekrna_dps[0][st][en][a]));
              goto loopend;
            }
          }
          if (!random_model && (a == DP_P || a == DP_U)) {
            energy_t rnastructureval = a == DP_P ?
                rnastructure_dp.v.f(st + 1, en + 1) : rnastructure_dp.w.f(st + 1, en + 1);
            if (((kekrna0 < constants::CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
                (kekrna0 < constants::CAP_E && kekrna0 != rnastructureval))) {
              errors.push_back(sfmt("rnastructure at %d %d %d: %d != %d",
                  st, en, a, rnastructureval, kekrna_dps[0][st][en][a]));
              goto loopend;
            }
          }
        }
      }
    }
    loopend:
    return errors;
  }

  error_t KekrnaComputeAndCheckState() {
    error_t errors;
    // Kekrna.
    std::vector<energy_t> kekrna_ctd_efns;
    std::vector<energy_t> kekrna_optimal_efns;
    for (const auto& options : kekrnas) {
      Context ctx(r, em, options);
      auto computed = ctx.Fold();
      kekrna_dps.emplace_back(std::move(gdp));
      // First compute with the CTDs that fold returned to check the energy.
      kekrna_ctd_efns.push_back(energy::ComputeEnergyWithCtds(computed, *em).energy);
      // Also check that the optimal CTD configuration has the same energy.
      // Note that it might not be the same, so we can't do an equality check.
      kekrna_optimal_efns.push_back(energy::ComputeEnergy(computed.s, *em).energy);
      kekrna_computeds.push_back(std::move(computed));
    }

    // Check kekrna energies.
    for (int i = 0; i < int(kekrnas.size()); ++i) {
      if (kekrna_computeds[0].energy != kekrna_computeds[i].energy ||
          kekrna_computeds[0].energy != kekrna_ctd_efns[i] ||
          kekrna_computeds[0].energy != kekrna_optimal_efns[i])
        errors.push_back(sfmt("kekrna %d: %d (dp) %d (ctd efn) %d (efn) != mfe %d",
            kekrna_computeds[i].energy, kekrna_ctd_efns[i], kekrna_optimal_efns[i], kekrna_computeds[0].energy));
    }

    return errors;
  }

  error_t RnastructureComputeAndCheckState() {
    error_t errors;
    auto rnastructure_computed = rnastructure.FoldAndDpTable(r, &rnastructure_dp);
    auto rnastructure_efn = rnastructure.Efn(rnastructure_computed.s);
    if (kekrna_computeds[0].energy != rnastructure_computed.energy ||
        kekrna_computeds[0].energy != rnastructure_efn)
      errors.push_back(sfmt("mfe: rnastructure %d (dp), %d (efn) != mfe %d",
          rnastructure_computed.energy, rnastructure_efn, kekrna_computeds[0].energy));
    return errors;
  }

  error_t CheckBruteForce() {
    error_t errors;
    context_options_t options(context_options_t::TableAlg::TWO, context_options_t::SuboptimalAlg::ZERO);
    Context ctx(r, em, options);
    auto subopt_brute = FoldBruteForce(r, *em, SUBOPT_BRUTE_MAX_STRUCTURES);
    auto subopt_kekrna = ctx.Suboptimal(-1, SUBOPT_BRUTE_MAX_STRUCTURES);

    AppendErrors(errors, MaybePrependHeader(CheckSuboptimal(subopt_brute, true), "brute suboptimal:"));
    AppendErrors(errors, MaybePrependHeader(CheckSuboptimal(subopt_kekrna, true), "kekrna suboptimal:"));
    AppendErrors(errors, MaybePrependHeader(CheckSuboptimalPair(subopt_brute, subopt_kekrna), "brute vs kekrna suboptimal:"));
    return errors;
  }

};


int main(int argc, char* argv[]) {
  std::mt19937 eng(uint_fast32_t(time(nullptr)));
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("-1")},
      {"random", ArgParse::option_t("use random energy models (disables comparison to RNAstructure)")}
  });
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 2, "require min and max length");
  const int min_len = atoi(pos[0].c_str());
  const int max_len = atoi(pos[1].c_str());
  verify_expr(min_len > 0, "invalid min length");
  verify_expr(max_len >= min_len, "invalid max len");

  bridge::Rnastructure rnastructure("extern/rnark/data_tables/", false);
  std::vector<context_options_t> kekrnas;
  for (auto table_alg : context_options_t::TABLE_ALGS)
    kekrnas.emplace_back(table_alg);

  auto start_time = std::chrono::steady_clock::now();
  const auto interval = atoi(argparse.GetOption("print-interval").c_str());
  std::uniform_int_distribution<int> len_dist(min_len, max_len);
  const auto t04em = energy::LoadEnergyModelFromArgParse(argparse);
  const bool random_model = argparse.HasFlag("random");
  for (int64_t i = 0;; ++i) {
    if (interval > 0 && std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start_time).count() > interval) {
      printf("Fuzzed %" PRId64" RNA\n", i);
      start_time = std::chrono::steady_clock::now();
    }
    int len = len_dist(eng);
    auto r = GenerateRandomPrimary(len, eng);

    uint_fast32_t seed = eng();
    Fuzzer fuzzer(r, t04em, random_model, seed, kekrnas, rnastructure);
    const auto res = fuzzer.Run();
    if (!res.empty()) {
      for (const auto& s : res)
        printf("%s\n", s.c_str());
      printf("\n");
    }
  }
}

