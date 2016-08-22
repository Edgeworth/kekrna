#include <cstdio>
#include <random>
#include <chrono>
#include "bridge/bridge.h"
#include "parsing.h"
#include "energy/energy_model.h"

using namespace kekrna;

rna_t GenerateRandomRna(int length) {
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < length; ++i)
    rna[i] = rand() % 4;
  return rna;
}

void FuzzRna(const rna_t& rna, bool use_random_energy_model,
    const std::vector<bridge::Kekrna>& kekrnas, const bridge::Rnastructure& rnastructure) {
  // Initialise everything here because we use goto later.
  int seed = rand();
  bool use_brute = rna.size() <= 22;
  bool dp_table_diff = false;
  int N = int(rna.size());
  dp_state_t rnastructure_state;
  folded_rna_t rnastructure_frna, brute_frna;
  energy_t rnastructure_efn = 0, brute_efn = 0;
  int st = 0, en = 0, a = 0;

  if (use_random_energy_model)
    energy::LoadRandomEnergyModel(seed);
  std::vector<fold::fold_state_t> kekrna_states;
  std::vector<folded_rna_t> kekrna_folds;
  std::vector<energy_t> kekrna_efns;
  for (const auto& kekrna : kekrnas) {
    fold::fold_state_t state;
    auto frna = kekrna.FoldAndDpTable(rna, &state);
    kekrna_efns.push_back(energy::ComputeEnergy(frna));
    kekrna_states.push_back(std::move(state));
    kekrna_folds.push_back(std::move(frna));
  }

  for (int i = 0; i < int(kekrnas.size()); ++i) {
    if (kekrna_folds[0].energy != kekrna_folds[i].energy ||
        kekrna_folds[0].energy != kekrna_efns[i])
      goto print_diff;
  }

  if (use_brute) {
    brute_frna = fold::FoldBruteForce(rna, nullptr);
    if (kekrna_folds[0].energy != brute_frna.energy)
      goto print_diff;
  }

  if (!use_random_energy_model) {
    rnastructure_frna = rnastructure.FoldAndDpTable(rna, &rnastructure_state);
    rnastructure_efn = rnastructure.Efn(rnastructure_frna);
    if (kekrna_folds[0].energy != rnastructure_frna.energy ||
        kekrna_folds[0].energy != rnastructure_efn)
      goto print_diff;
  }

  for (st = N - 1; st >= 0; --st) {
    for (en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (a = 0; a < fold::DP_SIZE; ++a) {
        auto kekrna0 = kekrna_states[0].dp_table[st][en][a];
        for (int i = 0; i < int(kekrnas.size()); ++i) {
          auto kekrnai = kekrna_states[i].dp_table[st][en][a];
          // If meant to be infinity and not.
          if (((kekrna0 < constants::CAP_E) != (kekrnai < constants::CAP_E)) ||
              (kekrna0 < constants::CAP_E && kekrna0 != kekrnai)) {
            dp_table_diff = true;
            goto print_diff;
          }
        }
        if (!use_random_energy_model) {
          energy_t rnastructureval = constants::MAX_E;
          if (a == fold::DP_P)
            rnastructureval = rnastructure_state.v.f(st + 1, en + 1);
          else if (a == fold::DP_U)
            rnastructureval = rnastructure_state.w.f(st + 1, en + 1);
          if (rnastructureval != constants::MAX_E &&
              ((kekrna0 < constants::CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
              (kekrna0 < constants::CAP_E && kekrna0 != rnastructureval))) {
            dp_table_diff = true;
            goto print_diff;
          }
        }
      }
    }
  }
  return;
  print_diff:
  printf("Difference on len %d RNA %s\n", int(rna.size()), parsing::RnaToString(rna).c_str());
  if (use_random_energy_model)
    printf("  Using random energy model with seed: %d\n", seed);
  else
    printf("  Using T04 energy model.\n");
  for (int i = 0; i < int(kekrnas.size()); ++i) {
    printf("  Fold%d: %d (dp), %d (efn) - %s\n", i, kekrna_folds[i].energy, kekrna_efns[i],
        parsing::PairsToDotBracket(kekrna_folds[i].p).c_str());
  }
  if (use_brute)
    printf("  BruteFold: %d (dp), %d (efn) - %s\n", brute_frna.energy, brute_efn,
        parsing::PairsToDotBracket(brute_frna.p).c_str());
  if (!use_random_energy_model)
    printf("  RNAstructure: %d (dp), %d (efn) - %s\n", rnastructure_frna.energy,
        rnastructure_efn, parsing::PairsToDotBracket(rnastructure_frna.p).c_str());
  if (dp_table_diff) {
    printf("  DP table difference at %d %d %d:\n", st, en, a);
    for (int i = 0; i < int(kekrnas.size()); ++i)
      printf("    Fold%d: %d\n", i, kekrna_states[i].dp_table[st][en][a]);
    if (!use_random_energy_model)
      printf("    RNAstructure: V: %d W: %d\n",
          rnastructure_state.v.f(st + 1, en + 1),
          rnastructure_state.w.f(st + 1, en + 1));
  }
}

int main(int argc, char* argv[]) {
  srand(static_cast<unsigned int>(time(NULL)));
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("-1")},
      {"random", ArgParse::option_t("use random energy models (disables comparison to RNAstructure)")}
  });
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(
      pos.size() == 2,
      "require size and variance");
  int base_len = atoi(pos[0].c_str());
  int variance = atoi(pos[1].c_str());
  verify_expr(base_len > 0, "invalid length");
  verify_expr(variance >= 0, "invalid variance");

  bridge::Rnastructure rnastructure("extern/rnark/data_tables/", false);
  std::vector<bridge::Kekrna> kekrnas;
  for (const auto& fold_fn : fold::FOLD_FUNCTIONS) {
    kekrnas.emplace_back(fold_fn);
  }

  auto start_time = std::chrono::steady_clock::now();
  auto interval = atoi(argparse.GetOption("print-interval").c_str());
  for (int i = 0;; ++i) {
    int length = base_len;
    if (variance) length += rand() % variance;
    if (interval > 0 && std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start_time).count() > interval) {
      printf("Fuzzed %d RNA\n", i);
      start_time = std::chrono::steady_clock::now();
    }
    auto rna = GenerateRandomRna(length);
    FuzzRna(rna, argparse.HasFlag("random"), kekrnas, rnastructure);
  }
}

