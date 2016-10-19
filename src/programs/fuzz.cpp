// Copyright 2016, Eliot Courtney.
//
// This file is part of kekrna.
//
// kekrna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// kekrna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with kekrna.
// If not, see <http://www.gnu.org/licenses/>.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <random>
#include <set>
#include "bridge/bridge.h"
#include "bridge/kekrna.h"
#include "bridge/rnastructure.h"
#include "energy/load_model.h"
#include "energy/structure.h"
#include "fold/brute_fold.h"
#include "parsing.h"

using namespace kekrna;
using namespace fold;
using namespace fold::internal;

const int SUBOPT_BRUTE_MAX_STRUCTURES = 10000;
const energy_t SUBOPT_MAX_DELTA = 6;  // Same as RNAstructure default.

class Fuzzer {
public:
  typedef std::deque<std::string> error_t;

  Fuzzer(primary_t r_, const energy::EnergyModelPtr em_, bool random_model_, uint_fast32_t seed_,
      const bridge::Rnastructure& rnastructure_, bool do_subopt_, bool do_subopt_rnastructure_,
      int max_structures_, int brute_cutoff_)
      : r(std::move(r_)), em(random_model_ ? energy::LoadRandomEnergyModel(seed_) : em_),
        random_model(random_model_), seed(seed_), rnastructure(rnastructure_),
        do_subopt(do_subopt_), do_subopt_rnastructure(do_subopt_rnastructure_),
        max_structures(max_structures_), brute_cutoff(brute_cutoff_) {
    verify_expr(!do_subopt_rnastructure || do_subopt,
        "suboptimal folding testing must be enabled to test rnastructure suboptimal folding");
    verify_expr(!(random_model && do_subopt_rnastructure),
        "cannot use a random energy model with rnastructure");
  }

  error_t Run() {
    error_t errors;
    AppendErrors(errors, MaybePrependHeader(KekrnaComputeAndCheckState(), "kekrna:"));
    if (!random_model)
      AppendErrors(errors, MaybePrependHeader(RnastructureComputeAndCheckState(), "rnastructure:"));
    if (int(r.size()) <= brute_cutoff)
      AppendErrors(errors, MaybePrependHeader(CheckBruteForce(), "brute force:"));
    AppendErrors(errors, MaybePrependHeader(CheckDpTables(), "dp tables:"));
    if (do_subopt) AppendErrors(errors, MaybePrependHeader(CheckSuboptimal(), "suboptimal:"));

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
  const primary_t r;
  const energy::EnergyModelPtr em;
  const bool random_model;
  const uint_fast32_t seed;
  const bridge::Rnastructure& rnastructure;
  const bool do_subopt;
  const bool do_subopt_rnastructure;
  const int max_structures;
  const int brute_cutoff;

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
      if (suboptimal_set.count(val)) return true;
      suboptimal_set.insert(val);
    }
    return false;
  }

  error_t CheckSuboptimalResult(const std::vector<computed_t>& subopt, bool has_ctds) {
    error_t errors;
    // Check at least one suboptimal structure.
    if (subopt.empty()) errors.push_back("no structures returned");
    // Check MFE.
    if (!subopt.empty() && kekrna_computeds[0].energy != subopt[0].energy)
      errors.push_back(sfmt(
          "lowest structure energy %d != mfe %d", subopt[0].energy, kekrna_computeds[0].energy));

    // Only ones with CTDs set can do these tests.
    if (has_ctds) {
      // Check for duplicate structures.
      if (HasDuplicates(subopt)) errors.push_back("has duplicates");

      for (int i = 0; i < int(subopt.size()); ++i) {
        const auto& structure = subopt[i];
        auto suboptimal_efn = energy::ComputeEnergyWithCtds(structure, *em);
        if (suboptimal_efn.energy != structure.energy) {
          errors.push_back(sfmt(
              "structure %d: energy %d != efn %d", i, structure.energy, suboptimal_efn.energy));
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

  error_t CheckSuboptimalResultPair(
      const std::vector<computed_t>& a, const std::vector<computed_t>& b) {
    error_t errors;
    if (a.size() != b.size()) {
      errors.push_back(sfmt(
          "first has %d structures != second has %d structures", int(a.size()), int(b.size())));
    } else {
      for (int i = 0; i < int(a.size()); ++i) {
        if (a[i].energy != b[i].energy) {
          errors.push_back(
              sfmt("structure %d: first %d != second %d", i, a[i].energy, b[i].energy));
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
      for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        for (a = 0; a < DP_SIZE; ++a) {
          const auto kekrna0 = kekrna_dps[0][st][en][a];
          for (int i = 0; i < int(kekrna_dps.size()); ++i) {
            const auto kekrnai = kekrna_dps[i][st][en][a];
            // If meant to be infinity and not.
            if (((kekrna0 < CAP_E) != (kekrnai < CAP_E)) ||
                (kekrna0 < CAP_E && kekrna0 != kekrnai)) {
              errors.push_back(sfmt("kekrna %d at %d %d %d: %d != %d", i, st, en, a,
                  kekrna_dps[i][st][en][a], kekrna_dps[0][st][en][a]));
              goto loopend;
            }
          }
          if (!random_model && (a == DP_P || a == DP_U)) {
            energy_t rnastructureval = a == DP_P ? rnastructure_dp.v.f(st + 1, en + 1)
                : rnastructure_dp.w.f(st + 1, en + 1);
            if (((kekrna0 < CAP_E) != (rnastructureval < INFINITE_ENERGY - 1000) ||
                (kekrna0 < CAP_E && kekrna0 != rnastructureval))) {
              errors.push_back(sfmt("rnastructure at %d %d %d: %d != %d", st, en, a,
                  rnastructureval, kekrna_dps[0][st][en][a]));
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
    for (auto table_alg : context_opt_t::TABLE_ALGS) {
      Context ctx(r, em, context_opt_t(table_alg));
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
    for (int i = 0; i < int(kekrna_dps.size()); ++i) {
      if (kekrna_computeds[0].energy != kekrna_computeds[i].energy ||
          kekrna_computeds[0].energy != kekrna_ctd_efns[i] ||
          kekrna_computeds[0].energy != kekrna_optimal_efns[i])
        errors.push_back(
            sfmt("kekrna %d: %d (dp) %d (ctd efn) %d (efn) != mfe %d", kekrna_computeds[i].energy,
                kekrna_ctd_efns[i], kekrna_optimal_efns[i], kekrna_computeds[0].energy));
    }

    return errors;
  }

  error_t RnastructureComputeAndCheckState() {
    error_t errors;
    auto rnastructure_computed = rnastructure.FoldAndDpTable(r, rnastructure_dp);
    auto rnastructure_efn = rnastructure.Efn(rnastructure_computed.s);
    if (kekrna_computeds[0].energy != rnastructure_computed.energy ||
        kekrna_computeds[0].energy != rnastructure_efn)
      errors.push_back(sfmt("mfe: rnastructure %d (dp), %d (efn) != mfe %d",
          rnastructure_computed.energy, rnastructure_efn, kekrna_computeds[0].energy));
    return errors;
  }

  error_t CheckBruteForce() {
    error_t errors;
    context_opt_t options(
        context_opt_t::TableAlg::TWO, context_opt_t::SuboptimalAlg::ZERO);
    Context ctx(r, em, options);
    auto subopt_brute = SuboptimalBruteForce(r, *em, SUBOPT_BRUTE_MAX_STRUCTURES);
    auto subopt_kekrna = ctx.SuboptimalIntoVector(true, -1, SUBOPT_BRUTE_MAX_STRUCTURES);

    AppendErrors(
        errors, MaybePrependHeader(CheckSuboptimalResult(subopt_brute, true), "brute suboptimal:"));
    AppendErrors(errors,
        MaybePrependHeader(CheckSuboptimalResult(subopt_kekrna, true), "kekrna suboptimal:"));
    AppendErrors(errors, MaybePrependHeader(CheckSuboptimalResultPair(subopt_brute, subopt_kekrna),
        "brute vs kekrna suboptimal:"));
    return errors;
  }

  error_t CheckSuboptimal() {
    error_t errors;
    std::vector<std::vector<computed_t>> kekrna_subopts_delta, kekrna_subopts_num;
    for (auto subopt_alg : context_opt_t::SUBOPTIMAL_ALGS) {
      context_opt_t options(context_opt_t::TableAlg::TWO, subopt_alg);
      Context ctx(r, em, options);
      kekrna_subopts_delta.push_back(ctx.SuboptimalIntoVector(true, SUBOPT_MAX_DELTA, -1));
      kekrna_subopts_num.push_back(ctx.SuboptimalIntoVector(true, -1, max_structures));
    }

    for (int i = 0; i < int(kekrna_subopts_delta.size()); ++i) {
      AppendErrors(errors, MaybePrependHeader(CheckSuboptimalResult(kekrna_subopts_delta[i], true),
          sfmt("kekrna delta suboptimal %d:", i)));
      AppendErrors(errors,
          MaybePrependHeader(
              CheckSuboptimalResultPair(kekrna_subopts_delta[0], kekrna_subopts_delta[i]),
              sfmt("kekrna 0 vs kekrna %d delta suboptimal:", i)));
    }

    for (int i = 0; i < int(kekrna_subopts_num.size()); ++i) {
      AppendErrors(errors, MaybePrependHeader(CheckSuboptimalResult(kekrna_subopts_num[i], true),
          sfmt("kekrna num suboptimal %d:", i)));
      AppendErrors(errors,
          MaybePrependHeader(
              CheckSuboptimalResultPair(kekrna_subopts_num[0], kekrna_subopts_num[i]),
              sfmt("kekrna 0 vs kekrna %d num suboptimal:", i)));
    }

    if (do_subopt_rnastructure) {
      // Suboptimal folding. Ignore ones with MFE >= -SUBOPT_MAX_DELTA because RNAstructure does
      // strange things
      // when the energy for suboptimal structures is 0 or above.
      if (kekrna_computeds[0].energy < -SUBOPT_MAX_DELTA) {
        context_opt_t options(
            context_opt_t::TableAlg::TWO, context_opt_t::SuboptimalAlg::ONE);
        Context ctx(r, em, options);
        auto kekrna_subopt = ctx.SuboptimalIntoVector(true, SUBOPT_MAX_DELTA, -1);
        const auto rnastructure_subopt = rnastructure.SuboptimalIntoVector(r, SUBOPT_MAX_DELTA);
        AppendErrors(errors,
            MaybePrependHeader(CheckSuboptimalResult(kekrna_subopt, true), "kekrna suboptimal:"));
        AppendErrors(errors, MaybePrependHeader(CheckSuboptimalResult(rnastructure_subopt, false),
            "rnastructure suboptimal:"));
        AppendErrors(errors,
            MaybePrependHeader(CheckSuboptimalResultPair(kekrna_subopt, rnastructure_subopt),
                "kekrna vs rnastructure suboptimal:"));
      }
    }

    return errors;
  }
};

int main(int argc, char* argv[]) {
  std::mt19937 eng(uint_fast32_t(time(nullptr)));
  ArgParse argparse({{"print-interval",
      opt_t("status update every n seconds").Arg("-1")},
      {"random", opt_t("use random energy models (disables comparison to RNAstructure)")},
      {"no-subopt", opt_t("do not test suboptimal folding")},
      {"subopt-rnastructure", opt_t("test rnastructure suboptimal folding")},
      {"subopt-max-structures", opt_t(
          "maximum number of substructures for max-delta fuzz").Arg("5000")},
      {"afl", opt_t("reads one rna from stdin and fuzzes - useful for use with afl")},
      {"brute-cutoff", opt_t("maximum rna size to run brute force on").Arg("25")}
  });
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);

  const bridge::Rnastructure rnastructure("extern/miles_rnastructure/data_tables/", false);
  const auto t04em = energy::LoadEnergyModelFromArgParse(argparse);
  const bool random_model = argparse.HasFlag("random");
  const bool do_subopt = !argparse.HasFlag("no-subopt");
  const bool do_subopt_rnastructure = argparse.HasFlag("subopt-rnastructure");
  const int max_structures = atoi(argparse.GetOption("subopt-max-structures").c_str());
  const bool afl_mode = argparse.HasFlag("afl");
  const int brute_cutoff = atoi(argparse.GetOption("brute-cutoff").c_str());

  if (afl_mode) {
// AFL mode.
#ifdef __AFL_HAVE_MANUAL_CONTROL
    __AFL_INIT();
    while (__AFL_LOOP(1000)) {
#endif
    std::string data;
    std::size_t len;
    char buf[4096];
    while ((len = fread(buf, 1, sizeof(buf), stdin)) > 0)
      data += std::string(buf, len);
    if (data.size() > 0) {
      uint_fast32_t seed = eng();
      // Disable brute force testing for AFL since it's too slow.
      Fuzzer fuzzer(parsing::StringToPrimary(data), t04em, random_model, seed,
          rnastructure, do_subopt, do_subopt_rnastructure, max_structures, 0);
      const auto res = fuzzer.Run();
      if (!res.empty())
        abort();
    }
#ifdef __AFL_HAVE_MANUAL_CONTROL
    }
#endif
  } else {
    auto pos = argparse.GetPositional();
    verify_expr(pos.size() == 2, "require min and max length");
    const int min_len = atoi(pos[0].c_str());
    const int max_len = atoi(pos[1].c_str());
    const auto interval = atoi(argparse.GetOption("print-interval").c_str());

    verify_expr(min_len > 0, "invalid min length");
    verify_expr(max_len >= min_len, "invalid max len");
    std::uniform_int_distribution<int> len_dist(min_len, max_len);

    printf("Fuzzing [%d, %d] len RNAs - ", min_len, max_len);
    if (random_model)
      printf("random energy models");
    else
      printf("T04 energy model");
    if (do_subopt) printf(" - testing suboptimal folders");
    if (do_subopt_rnastructure) printf(" (including rnastructure)");
    printf("\n");

    // Normal mode.
    auto start_time = std::chrono::steady_clock::now();
    for (int64_t i = 0;; ++i) {
      if (interval > 0 &&
          std::chrono::duration_cast<std::chrono::seconds>(
              std::chrono::steady_clock::now() - start_time)
              .count() > interval) {
        printf("Fuzzed %" PRId64 " RNA\n", i);
        start_time = std::chrono::steady_clock::now();
      }
      int len = len_dist(eng);
      auto r = GenerateRandomPrimary(len, eng);

      uint_fast32_t seed = eng();
      Fuzzer fuzzer(r, t04em, random_model, seed, rnastructure, do_subopt, do_subopt_rnastructure,
          max_structures, brute_cutoff);
      const auto res = fuzzer.Run();
      if (!res.empty()) {
        for (const auto& s : res)
          printf("%s\n", s.c_str());
        printf("\n");
      }
    }
  }
}
