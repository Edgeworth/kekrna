#include "parsing.h"

namespace kekrna {
namespace parsing {

rna_t ParseRnaFromString(const std::string& s) {
  rna_t rna(s.size());
  for (int i = 0; i < s.size(); ++i) {
    rna[i] = CharToBase(s[i]);
  }
  return rna;
}

folded_rna_t ParseViennaRna(const std::string& rna_str, const std::string& pairs_str) {
  rna_t rna = ParseRnaFromString(rna_str);
  std::vector<int> pairs(rna_str.size(), -1);
  std::stack<int> s;
  for (int i = 0; i < rna_str.size(); ++i) {
    if (pairs_str[i] == '(') {
      s.push(i);
    } else if (pairs_str[i] == ')') {
      assert(!s.empty());
      pairs[i] = s.top();
      pairs[s.top()] = i;
      s.pop();
    }
  }
  return {rna, pairs};
}


void Parse2x2FromFile(const std::string& filename, energy::energy_t (&output)[4][4][4][4]) {
  FILE* fp = fopen(filename.c_str(), "r");
  for (int i = 0; i < 256; ++i) {
    base_t a = CharToBase((char) fgetc(fp));
    base_t b = CharToBase((char) fgetc(fp));
    base_t c = CharToBase((char) fgetc(fp));
    base_t d = CharToBase((char) fgetc(fp));
    assert(a != -1 && b != -1 && c != -1 && d != -1);
    int res = fscanf(fp, " %d ", &output[a][b][c][d]);
    assert(res == 1);
  }
  fclose(fp);
}

void ParseMapFromFile(const std::string& filename, std::unordered_map<std::string, energy::energy_t>& output) {
  FILE* fp = fopen(filename.c_str(), "r");
  char buf[1024];
  energy::energy_t energy;
  while (fscanf(fp, " %s %d ", buf, &energy) == 2)
    output[buf] = energy;
  fclose(fp);
}

void ParseInitiationEnergyFromFile(const std::string& filename, energy::energy_t (&output)[INITIATION_CACHE_SZ]) {
  FILE* fp = fopen(filename.c_str(), "r");
  energy::energy_t energy;
  int idx;
  while (fscanf(fp, "%d %d ", &idx, &energy) == 2)
    output[idx] = energy;
  fclose(fp);
}

void ParseHairpinMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(
      fp, "%d %d %d %d %d %d",
      &hairpin_uu_ga_first_mismatch, &hairpin_gg_first_mismatch,
      &hairpin_special_gu_closure, &hairpin_c3_loop, &hairpin_all_c_a, &hairpin_all_c_b);
  assert(res == 6);
  fclose(fp);
}

void ParseBulgeMiscDataFromFile(const std::string& filename) {
  FILE* fp = fopen(filename.c_str(), "r");
  int res = fscanf(fp, "%d", &bulge_special_c);
  assert(res == 1);
  fclose(fp);
}


}
}
