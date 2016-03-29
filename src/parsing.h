#ifndef KEKRNA_PARSING_H
#define KEKRNA_PARSING_H

#include <string>
#include <stack>
#include "base.h"

namespace kekrna {
namespace parsing {

rna_t ParseRnaFromString(const std::string& s);

folded_rna_t ParseViennaRna(const std::string& rna_str, const std::string& pairs_str);

void ParseStackingEnergiesFromFile(const std::string& filename);

}
}

#endif //KEKRNA_PARSING_H
