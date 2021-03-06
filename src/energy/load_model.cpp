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
#include "load_model.h"
#include "parsing.h"

namespace kekrna {
namespace energy {

namespace {
const energy_t RAND_MIN_ENERGY = -1000;
const energy_t RAND_MAX_ENERGY = 1000;
const int RAND_MAX_HAIRPIN_SZ = 8;
const int RAND_MAX_NUM_HAIRPIN = 50;

void Parse2x2FromFile(const std::string& filename, energy_t (& output)[4][4][4][4]) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    const base_t a = CharToBase(char(fgetc(fp)));
    const base_t b = CharToBase(char(fgetc(fp)));
    const base_t c = CharToBase(char(fgetc(fp)));
    const base_t d = CharToBase(char(fgetc(fp)));
    if (a == -1) break;
    verify_expr(a != -1 && b != -1 && c != -1 && d != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &output[a][b][c][d]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseMapFromFile(
    const std::string& filename, std::unordered_map<std::string, energy_t>& output) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  char buf[1024];
  energy_t energy;
  while (fscanf(fp, " %s %d ", buf, &energy) == 2)
    output[buf] = energy;
  fclose(fp);
}

void ParseInitiationEnergyFromFile(
    const std::string& filename, energy_t (& output)[energy::EnergyModel::INITIATION_CACHE_SZ]) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  energy_t energy;
  int idx;
  while (fscanf(fp, "%d %d ", &idx, &energy) == 2) {
    verify_expr(idx < energy::EnergyModel::INITIATION_CACHE_SZ,
        "out of bounds index in %s", filename.c_str());
    output[idx] = energy;
  }
  fclose(fp);
}

void ParseInternalLoop1x1FromFile(const std::string& filename, energy::EnergyModel& em) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    const base_t a = CharToBase(char(fgetc(fp)));
    const base_t b = CharToBase(char(fgetc(fp)));
    const base_t c = CharToBase(char(fgetc(fp)));
    const base_t d = CharToBase(char(fgetc(fp)));
    const base_t e = CharToBase(char(fgetc(fp)));
    const base_t f = CharToBase(char(fgetc(fp)));
    if (a == -1) break;
    verify_expr(a != -1 && b != -1 && c != -1 && d != -1 && e != -1 && f != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &em.internal_1x1[a][b][c][d][e][f]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseInternalLoop1x2FromFile(const std::string& filename, energy::EnergyModel& em) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    const base_t a = CharToBase(char(fgetc(fp)));
    const base_t b = CharToBase(char(fgetc(fp)));
    const base_t c = CharToBase(char(fgetc(fp)));
    const base_t d = CharToBase(char(fgetc(fp)));
    const base_t e = CharToBase(char(fgetc(fp)));
    const base_t f = CharToBase(char(fgetc(fp)));
    const base_t g = CharToBase(char(fgetc(fp)));
    if (a == -1) break;
    verify_expr(
        a != -1 && b != -1 && c != -1 && d != -1 && e != -1 && f != -1 && g != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &em.internal_1x2[a][b][c][d][e][f][g]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseInternalLoop2x2FromFile(const std::string& filename, energy::EnergyModel& em) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    const base_t a = CharToBase(char(fgetc(fp)));
    const base_t b = CharToBase(char(fgetc(fp)));
    const base_t c = CharToBase(char(fgetc(fp)));
    const base_t d = CharToBase(char(fgetc(fp)));
    const base_t e = CharToBase(char(fgetc(fp)));
    const base_t f = CharToBase(char(fgetc(fp)));
    const base_t g = CharToBase(char(fgetc(fp)));
    const base_t h = CharToBase(char(fgetc(fp)));
    if (a == -1) break;
    verify_expr(
        a != -1 && b != -1 && c != -1 && d != -1 && e != -1 && f != -1 && g != -1 && h != -1,
        "expected base");
    verify_expr(
        fscanf(fp, " %d ", &em.internal_2x2[a][b][c][d][e][f][g][h]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseDangleDataFromFile(const std::string& filename, energy_t (& output)[4][4][4]) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");
  while (1) {
    base_t a = CharToBase(char(fgetc(fp)));
    base_t b = CharToBase(char(fgetc(fp)));
    base_t c = CharToBase(char(fgetc(fp)));
    if (a == -1) break;
    verify_expr(a != -1 && b != -1 && c != -1, "expected base");
    verify_expr(fscanf(fp, " %d ", &output[a][b][c]) == 1, "expected energy");
  }
  fclose(fp);
}

void ParseMiscDataFromFile(const std::string& filename, energy::EnergyModel& em) {
  FILE* fp = fopen(filename.c_str(), "r");
  verify_expr(fp != nullptr, "could not open file");

#define READ_DATA(var)                                                      \
  do {                                                                      \
    while (1) {                                                             \
      std::string line = sgetline(fp);                                      \
      verify_expr(line.size() > 0, "unexpected EOF or error");              \
      if (line[0] == '/' || line[0] == '\n') continue;                      \
      verify_expr(sscanf(line.c_str(), "%d", &(var)) == 1, "expected int"); \
      break;                                                                \
    }                                                                       \
  } while (0)

  // Bulge loops.
  READ_DATA(em.bulge_special_c);

  // Coaxial stacking.
  READ_DATA(em.coax_mismatch_non_contiguous);
  READ_DATA(em.coax_mismatch_wc_bonus);
  READ_DATA(em.coax_mismatch_gu_bonus);

  // Hairpin loops.
  READ_DATA(em.hairpin_uu_ga_first_mismatch);
  READ_DATA(em.hairpin_gg_first_mismatch);
  READ_DATA(em.hairpin_special_gu_closure);
  READ_DATA(em.hairpin_c3_loop);
  READ_DATA(em.hairpin_all_c_a);
  READ_DATA(em.hairpin_all_c_b);

  // Internal loops.
  READ_DATA(em.internal_asym);
  READ_DATA(em.internal_augu_penalty);

  // Multiloop data.
  READ_DATA(em.multiloop_hack_a);
  READ_DATA(em.multiloop_hack_b);

  // AU/GU penalty
  READ_DATA(em.augu_penalty);
#undef READ_DATA

  fclose(fp);
}
}

EnergyModelPtr LoadRandomEnergyModel(uint_fast32_t seed) {
  auto em = std::make_shared<EnergyModel>();
  std::mt19937 eng(seed);
  std::uniform_int_distribution<energy_t> energy_dist(RAND_MIN_ENERGY, RAND_MAX_ENERGY);
  std::uniform_int_distribution<energy_t> nonneg_energy_dist(0, RAND_MAX_ENERGY);
#define RANDOMISE_DATA(d)                                          \
  do {                                                             \
    auto dp = reinterpret_cast<energy_t*>(&(d));                   \
    for (unsigned int i = 0; i < sizeof(d) / sizeof(*dp); ++i) {   \
      dp[i] = energy_dist(eng);                                    \
    }                                                              \
  } while (0)

  RANDOMISE_DATA(em->stack);
  RANDOMISE_DATA(em->terminal);
  RANDOMISE_DATA(em->internal_init);
  RANDOMISE_DATA(em->internal_1x1);
  RANDOMISE_DATA(em->internal_1x2);
  RANDOMISE_DATA(em->internal_2x2);
  RANDOMISE_DATA(em->internal_2x3_mismatch);
  RANDOMISE_DATA(em->internal_other_mismatch);
  em->internal_asym =
      nonneg_energy_dist(eng);  // This needs to be non-negative for some optimisations.
  RANDOMISE_DATA(em->internal_augu_penalty);
  RANDOMISE_DATA(em->bulge_init);
  RANDOMISE_DATA(em->bulge_special_c);
  RANDOMISE_DATA(em->hairpin_init);
  RANDOMISE_DATA(em->hairpin_uu_ga_first_mismatch);
  RANDOMISE_DATA(em->hairpin_gg_first_mismatch);
  RANDOMISE_DATA(em->hairpin_special_gu_closure);
  RANDOMISE_DATA(em->hairpin_c3_loop);
  RANDOMISE_DATA(em->hairpin_all_c_a);
  RANDOMISE_DATA(em->hairpin_all_c_b);

  em->hairpin.clear();
  std::uniform_int_distribution<int> hairpin_size_dist(HAIRPIN_MIN_SZ, RAND_MAX_HAIRPIN_SZ);
  static_assert(HAIRPIN_MIN_SZ <= RAND_MAX_HAIRPIN_SZ,
      "HAIRPIN_MIN_SZ > RAND_MAX_HAIRPIN does not make sense");
  std::uniform_int_distribution<int> num_hairpin_dist(0, RAND_MAX_NUM_HAIRPIN);
  int num_hairpin = num_hairpin_dist(eng);
  for (int i = 0; i < num_hairpin; ++i) {
    auto hairpin = parsing::PrimaryToString(GenerateRandomPrimary(hairpin_size_dist(eng), eng));
    em->hairpin[hairpin] = energy_dist(eng);
  }

  RANDOMISE_DATA(em->multiloop_hack_a);
  RANDOMISE_DATA(em->multiloop_hack_b);
  RANDOMISE_DATA(em->dangle5);
  RANDOMISE_DATA(em->dangle3);
  RANDOMISE_DATA(em->coax_mismatch_non_contiguous);
  RANDOMISE_DATA(em->coax_mismatch_wc_bonus);
  RANDOMISE_DATA(em->coax_mismatch_gu_bonus);
  RANDOMISE_DATA(em->augu_penalty);
#undef RANDOMISE_DATA

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Correct things to be 180 degree rotations if required.
          em->stack[c][d][a][b] = em->stack[a][b][c][d];
          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              em->internal_1x1[d][e][f][a][b][c] = em->internal_1x1[a][b][c][d][e][f];
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  em->internal_2x2[e][f][g][h][a][b][c][d] =
                      em->internal_2x2[a][b][c][d][e][f][g][h];
                }
              }
            }
          }
        }
      }
    }
  }

  std::string reason;
  verify_expr(em->IsValid(&reason), "invalid energy model: %s", reason.c_str());
  return em;
}

EnergyModelPtr LoadEnergyModelFromDataDir(const std::string& data_dir) {
  auto em = std::make_shared<EnergyModel>();
  // Stacking interaction data.
  Parse2x2FromFile(data_dir + "/stacking.data", em->stack);

  // Terminal mismatch data.
  Parse2x2FromFile(data_dir + "/terminal.data", em->terminal);

  // Hairpin data.
  ParseMapFromFile(data_dir + "/hairpin.data", em->hairpin);
  ParseInitiationEnergyFromFile(data_dir + "/hairpin_initiation.data", em->hairpin_init);

  // Bulge loop data.
  ParseInitiationEnergyFromFile(data_dir + "/bulge_initiation.data", em->bulge_init);

  // Internal loop data.
  ParseInitiationEnergyFromFile(data_dir + "/internal_initiation.data", em->internal_init);
  ParseInternalLoop1x1FromFile(data_dir + "/internal_1x1.data", *em);
  ParseInternalLoop1x2FromFile(data_dir + "/internal_1x2.data", *em);
  ParseInternalLoop2x2FromFile(data_dir + "/internal_2x2.data", *em);
  Parse2x2FromFile(data_dir + "/internal_2x3_mismatch.data", em->internal_2x3_mismatch);
  Parse2x2FromFile(data_dir + "/internal_other_mismatch.data", em->internal_other_mismatch);

  // Dangle data.
  ParseDangleDataFromFile(data_dir + "/dangle3.data", em->dangle3);
  ParseDangleDataFromFile(data_dir + "/dangle5.data", em->dangle5);

  // Other misc data.
  ParseMiscDataFromFile(data_dir + "/misc.data", *em);

  std::string reason;
  verify_expr(em->IsValid(&reason), "invalid energy model: %s", reason.c_str());
  return em;
}

EnergyModelPtr LoadEnergyModelFromArgParse(const ArgParse& argparse) {
  if (argparse.HasFlag("seed")) {
    return LoadRandomEnergyModel(uint_fast32_t(atoi(argparse.GetOption("seed").c_str())));
  } else {
    return LoadEnergyModelFromDataDir(argparse.GetOption("data-path"));
  }
}
}
}
