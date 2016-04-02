#include "base.h"
#include "parsing.h"
#include "globals.h"

namespace kekrna {

void Init() {
  parsing::Parse2x2FromFile("data/stacking.data", stacking_e);
  parsing::Parse2x2FromFile("data/terminal.data", terminal_e);
  parsing::ParseMapFromFile("data/hairpin.data", hairpin_e);
  parsing::ParseInitiationEnergyFromFile("data/internal_initiation.data", internal_init);
  parsing::ParseInitiationEnergyFromFile("data/bulge_initiation.data", bulge_init);
  parsing::ParseInitiationEnergyFromFile("data/hairpin_initiation.data", hairpin_init);
  parsing::ParseHairpinMiscDataFromFile("data/hairpin_misc.data");
  parsing::ParseBulgeMiscDataFromFile("data/bulge_misc.data");
}

base_t CharToBase(char c) {
  switch (c) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'U':
      return U;
    default:
      return -1;
  }
}

char BaseToChar(base_t b) {
  if (b < 0 || b > 4) return '?';
  return "ACGU"[b];
}

}
