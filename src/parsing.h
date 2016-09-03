#ifndef KEKRNA_PARSING_H
#define KEKRNA_PARSING_H

#include <string>
#include <stack>
#include <unordered_map>
#include "base.h"
#include "energy/energy_globals.h"

namespace kekrna {
namespace parsing {

primary_t StringToPrimary(const std::string& s);
std::string PrimaryToString(const primary_t& r);
secondary_t ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str);
std::vector<int> DotBracketToPairs(const std::string& pairs_str);
std::string PairsToDotBracket(const std::vector<int>& pairs);
std::string ComputedToCtdString(const computed_t& computed);

void Parse2x2FromFile(const std::string& filename, energy_t (& output)[4][4][4][4]);
void ParseMapFromFile(const std::string& filename, std::unordered_map<std::string, energy_t>& output);
void ParseInitiationEnergyFromFile(const std::string& filename, energy_t (& output)[energy::INITIATION_CACHE_SZ]);
void ParseInternalLoop1x1FromFile(const std::string& filename);
void ParseInternalLoop1x2FromFile(const std::string& filename);
void ParseInternalLoop2x2FromFile(const std::string& filename);
void ParseDangleDataFromFile(const std::string& filename, energy_t (& output)[4][4][4]);
void ParseMiscDataFromFile(const std::string& filename);

}
}

#endif //KEKRNA_PARSING_H
