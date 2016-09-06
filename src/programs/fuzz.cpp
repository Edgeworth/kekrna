#include <cstdio>
#include <random>
#include <chrono>
#include <cinttypes>
#include <set>
#include "fold/brute_fold.h"
#include "bridge/bridge.h"
#include "parsing.h"
#include "energy/load_model.h"

using namespace kekrna;
using namespace fold;

const int SUBOPT_MAX_STRUCTURES = 50;

template<typename RandomEngine>
void FuzzRna(const primary_t& r, bool use_random_energy_model, const energy::EnergyModel& loaded,
    const std::vector<context_options_t>& kekrnas, const bridge::Rnastructure& rnastructure,
    RandomEngine& eng) {
  uint_fast32_t seed = eng();
  auto em = loaded;
  if (use_random_energy_model)
    em = energy::LoadRandomEnergyModel(seed);

  // Kekrna.
  std::vector<Context> kekrna_ctxs;
  std::vector<computed_t> kekrna_folds;
  std::vector<energy_t> kekrna_efns;
  std::vector<energy_t> kekrna_optimal_efns;
  for (const auto& options : kekrnas) {
    Context ctx(r, em, options);
    auto computed = ctx.Fold();
    // First compute with the CTDs that fold returned to check the energy.
    kekrna_efns.push_back(energy::ComputeEnergyWithCtds(computed, em).energy);
    // Also check that the optimal CTD configuration has the same energy.
    // Note that it might not be the same, so we can't do an equality check.
    kekrna_optimal_efns.push_back(energy::ComputeEnergy(computed.s, em).energy);
    kekrna_ctxs.push_back(std::move(ctx));
    kekrna_folds.push_back(std::move(computed));
  }

  // Check kekrna energies.
  bool mfe_diff = false;
  for (int i = 0; i < int(kekrnas.size()); ++i) {
    if (kekrna_folds[0].energy != kekrna_folds[i].energy ||
        kekrna_folds[0].energy != kekrna_efns[i] ||
        kekrna_folds[0].energy != kekrna_optimal_efns[i])
      mfe_diff = true;
  }

  // Brute force.
  bool use_brute = r.size() <= 24;
  computed_t brute_computed;
  if (use_brute) {
    brute_computed = FoldBruteForce(r, em);
    if (kekrna_folds[0].energy != brute_computed.energy)
      mfe_diff = true;
  }

  // RNAstructure.
  computed_t rnastructure_computed;
  energy_t rnastructure_efn = 0;
  dp_state_t rnastructure_state;
  if (!use_random_energy_model) {
    rnastructure_computed = rnastructure.FoldAndDpTable(r, &rnastructure_state);
    rnastructure_efn = rnastructure.Efn(rnastructure_computed.s);
    if (kekrna_folds[0].energy != rnastructure_computed.energy ||
        kekrna_folds[0].energy != rnastructure_efn)
      mfe_diff = true;
  }

  // Suboptimal folding:
  bool suboptimal_mfe_diff = false;
  context_options_t options(context_options_t::TableAlg::TWO,
      context_options_t::SuboptimalAlg::ZERO, -1, SUBOPT_MAX_STRUCTURES);
  Context ctx(r, em, options);
  auto computeds = ctx.Suboptimal();
  // Check MFE.
  if (kekrna_folds[0].energy != computeds[0].energy)
    suboptimal_mfe_diff = true;
  // Check for duplicate structures.
  bool suboptimal_duplicate = false;
  // Check efn gives the same value.
  bool suboptimal_efn_diff = false;
  // If energies are different but everything else is the same, it is still a bug.
  std::set<std::pair<secondary_t, std::vector<Ctd>>> suboptimal_set;
  for (const auto& computed : computeds) {
    auto suboptimal_efn = energy::ComputeEnergyWithCtds(computed, em);
    if (suboptimal_efn.energy != computed.energy)
      suboptimal_efn_diff = true;
    auto val = std::make_pair(computed.s, computed.base_ctds);
    if (suboptimal_set.count(val))
      suboptimal_duplicate = true;
    suboptimal_set.insert(std::move(val));
  }

  int st = 0, en = 0, a = 0;
  bool dp_table_diff = false;
  const int N = int(r.size());
  for (st = N - 1; st >= 0; --st) {
    for (en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (a = 0; a < DP_SIZE; ++a) {
        const auto kekrna0 = kekrna_ctxs[0].GetDpState()[st][en][a];
        for (int i = 0; i < int(kekrnas.size()); ++i) {
          const auto kekrnai = kekrna_ctxs[i].GetDpState()[st][en][a];
          // If meant to be infinity and not.
          if (((kekrna0 < constants::CAP_E) != (kekrnai < constants::CAP_E)) ||
              (kekrna0 < constants::CAP_E && kekrna0 != kekrnai)) {
            dp_table_diff = true;
            goto loop_end;
          }
        }
        if (!use_random_energy_model) {
          energy_t rnastructureval = constants::MAX_E;
          if (a == DP_P)
            rnastructureval = rnastructure_state.v.f(st + 1, en + 1);
          else if (a == DP_U)
            rnastructureval = rnastructure_state.w.f(st + 1, en + 1);
          if (rnastructureval != constants::MAX_E &&
              ((kekrna0 < constants::CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
                  (kekrna0 < constants::CAP_E && kekrna0 != rnastructureval))) {
            dp_table_diff = true;
            goto loop_end;
          }
        }
      }
    }
  }
  loop_end:;
  if (mfe_diff || dp_table_diff) {
    printf("Difference on len %zu RNA %s\n", r.size(), parsing::PrimaryToString(r).c_str());
    if (use_random_energy_model)
      printf("  Using random energy model with seed: %" PRIuFAST32 "\n", seed);
    else
      printf("  Using T04 energy model.\n");
    if (mfe_diff) {
      for (int i = 0; i < int(kekrnas.size()); ++i) {
        printf("  Fold%d: %d (dp), %d (efn), %d (optimal efn) - %s\n", i, kekrna_folds[i].energy, kekrna_efns[i],
            kekrna_optimal_efns[i], parsing::PairsToDotBracket(kekrna_folds[i].s.p).c_str());
      }
      if (use_brute)
        printf("  BruteFold: %d - %s\n", brute_computed.energy,
            parsing::PairsToDotBracket(brute_computed.s.p).c_str());
      if (!use_random_energy_model)
        printf("  RNAstructure: %d (dp), %d (optimal efn) - %s\n", rnastructure_computed.energy,
            rnastructure_efn, parsing::PairsToDotBracket(rnastructure_computed.s.p).c_str());
    }
    if (dp_table_diff) {
      printf("  DP table difference at %d %d %d:\n", st, en, a);
      for (int i = 0; i < int(kekrnas.size()); ++i)
        printf("    Fold%d: %d\n", i, kekrna_ctxs[i].GetDpState()[st][en][a]);
      if (!use_random_energy_model)
        printf("    RNAstructure: V: %d W: %d\n",
            rnastructure_state.v.f(st + 1, en + 1),
            rnastructure_state.w.f(st + 1, en + 1));
    }
    if (suboptimal_mfe_diff) {
      printf("  Suboptimal: First structure was not MFE structure.\n");
    }
    if (suboptimal_efn_diff) {
      printf("  Suboptimal: Energy differs to EFN.\n");
    }
    if (suboptimal_duplicate) {
      printf("  Suboptimal: Duplicate structure.\n");
    }
  }
}

int main(int argc, char* argv[]) {
  std::mt19937 eng(uint_fast32_t(time(nullptr)));
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("-1")},
      {"random", ArgParse::option_t("use random energy models (disables comparison to RNAstructure)")}
  });
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(
      pos.size() == 2,
      "require min and max length");
  int min_len = atoi(pos[0].c_str());
  int max_len = atoi(pos[1].c_str());
  verify_expr(min_len > 0, "invalid min length");
  verify_expr(max_len >= min_len, "invalid max len");

  bridge::Rnastructure rnastructure("extern/rnark/data_tables/", false);
  std::vector<context_options_t> kekrnas;
  for (auto table_alg : context_options_t::TABLE_ALGS) {
    kekrnas.emplace_back(table_alg);
  }

  auto start_time = std::chrono::steady_clock::now();
  auto interval = atoi(argparse.GetOption("print-interval").c_str());
  std::uniform_int_distribution<int> len_dist(min_len, max_len);
  auto em = energy::LoadEnergyModelFromArgParse(argparse);
  for (int64_t i = 0;; ++i) {
    int length = len_dist(eng);
    if (interval > 0 && std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start_time).count() > interval) {
      printf("Fuzzed %d RNA\n", i);
      start_time = std::chrono::steady_clock::now();
    }
    auto r = GenerateRandomPrimary(length, eng);
    FuzzRna(r, argparse.HasFlag("random"), em, kekrnas, rnastructure, eng);
  }
}

