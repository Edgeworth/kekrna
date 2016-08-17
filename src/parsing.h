#ifndef KEKRNA_PARSING_H
#define KEKRNA_PARSING_H

#include <string>
#include <stack>
#include <unordered_map>
#include "base.h"
#include "globals.h"

namespace kekrna {
namespace parsing {

rna_t StringToRna(const std::string& s);

std::string RnaToString(const rna_t& rna);

folded_rna_t ParseDotBracketRna(const std::string& rna_str, const std::string& pairs_str);

std::vector<int> DotBracketToPairs(const std::string& pairs_str);

std::string PairsToDotBracket(const std::vector<int>& pairs);


void Parse2x2FromFile(const std::string& filename, energy_t (& output)[4][4][4][4]);

void ParseMapFromFile(const std::string& filename, std::unordered_map<std::string, energy_t>& output);

void ParseInitiationEnergyFromFile(const std::string& filename, energy_t (& output)[INITIATION_CACHE_SZ]);

void ParseInternalLoop1x1FromFile(const std::string& filename);

void ParseInternalLoop1x2FromFile(const std::string& filename);

void ParseInternalLoop2x2FromFile(const std::string& filename);

void ParseDangleDataFromFile(const std::string& filename, energy_t (& output)[4][4][4]);

void ParseMiscDataFromFile(const std::string& filename);

}
}

#endif //KEKRNA_PARSING_H
