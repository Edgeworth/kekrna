#include <cstdio>
#include <cassert>
#include <random>
#include <chrono>
#include "bridge/bridge.h"
#include "parsing.h"

using namespace kekrna;

void FuzzRnaOfLength(int length) {
  assert(length * 2 < 31);
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < 1 << (length * 2); ++i) {
    int data = i;
    for (int j = 0; j < length; ++j) {
      rna[j] = base_t(data & 0b11);
      data >>= 2;
    }
    energy_t a = fold::Fold0(rna).energy;
    energy_t a_efn = energy::ComputeEnergy();
    std::string a_db = parsing::PairsToDotBracket(p).c_str();
    energy_t b = fold::FoldBruteForce(rna).energy;
    energy_t b_efn = energy::ComputeEnergy();
    std::string b_db = parsing::PairsToDotBracket(p).c_str();

    if (a != b || a != a_efn || b != b_efn) {
      printf(
          "Diff. DP: %d (efn %d) != BF: %d (efn %d), on: %s\n%s vs %s\n",
          a, a_efn, b, b_efn, parsing::RnaToString(rna).c_str(), a_db.c_str(), b_db.c_str());
      return;
    }
  }
}

rna_t GenerateRandomRna(int length) {
  rna_t rna(std::size_t(length), 0);
  for (int i = 0; i < length; ++i)
    rna[i] = rand() % 4;
  return rna;
}

void FuzzRna(const rna_t& rna,
    const bridge::Kekrna& kekrna, const bridge::Rnastructure& rnastructure) {
  auto kekrna_dp = kekrna.Fold(rna);
  auto kekrna_efn = kekrna.Efn(kekrna_dp);
  auto rnastructure_dp = rnastructure.Fold(rna);
  auto rnastructure_efn = rnastructure.Efn(rnastructure_dp);

  if (kekrna_dp.energy != kekrna_efn ||
      kekrna_dp.energy != rnastructure_dp.energy ||
      kekrna_dp.energy != rnastructure_efn) {
    printf(
        "Diff on %s\n  Rnastructure: dp %d, efn %d\n    %s\n  Kekrna: dp %d, efn %d\n    %s\n\n",
        parsing::RnaToString(rna).c_str(), rnastructure_dp.energy, rnastructure_efn,
        parsing::PairsToDotBracket(rnastructure_dp.p).c_str(),
        kekrna_dp.energy, kekrna_efn, parsing::PairsToDotBracket(kekrna_dp.p).c_str()
    );
    return;
  }
}

void FuzzComputeTables(const rna_t& rna) {
  SetRna(rna);
  auto table0 = fold::ComputeTables0();
  auto table1 = fold::ComputeTables1();
  auto table2 = fold::ComputeTables2();
  auto table3 = fold::ComputeTables3();
  int N = int(rna.size());
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      for (int a = 0; a < fold::DP_SIZE; ++a) {
        if ((table0[st][en][a] != table1[st][en][a] ||
            table0[st][en][a] != table2[st][en][a] ||
            table0[st][en][a] != table3[st][en][a]) &&
            table0[st][en][a] < constants::CAP_E) {
          printf("Diff on %s\n %d %d %d: %d %d %d %d\n",
              parsing::RnaToString(rna).c_str(), st, en, a,
              table0[st][en][a], table1[st][en][a], table2[st][en][a],
              table3[st][en][a]);
          return;
        }
      }
    }
  }
}

// test types:
// against efn / brute fold, against rnastructure
// if against itself, try random energy models
// input types: stdin, exhaustive rnas (up to size x), random rnas (of size, and variance)
// usage: -e?, -f --stdin --brute --rand
// if brute force, up to size x
// if random, size, variance
// for efn, try all foldings as well?

int main(int argc, char* argv[]) {
  srand(static_cast<unsigned int>(time(NULL)));
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("60")}
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

  // bridge::Rnastructure rnastructure("extern/rnark/data_tables/", false);

  auto start_time = std::chrono::steady_clock::now();
  auto interval = atoi(argparse.GetOption("print-interval").c_str());
  for (int i = 0;; ++i) {
    int length = base_len;
    if (variance) length += rand() % variance;
    if (std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - start_time).count() > interval) {
      printf("Fuzzed %d RNA\n", i);
      start_time = std::chrono::steady_clock::now();
    }
    auto rna = GenerateRandomRna(length);
    FuzzComputeTables(rna);
  }
}

