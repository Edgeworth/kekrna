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
#include "structure.h"

namespace kekrna {
namespace energy {

const char* CtdToName(Ctd ctd) {
  switch (ctd) {
    case CTD_NA:
      return "n/a";
    case CTD_UNUSED:
      return "unused";
    case CTD_3_DANGLE:
      return "3' dangle";
    case CTD_5_DANGLE:
      return "5' dangle";
    case CTD_MISMATCH:
      return "terminal mismatch";
    case CTD_LCOAX_WITH_NEXT:
      return "left mismatch coax with next";
    case CTD_LCOAX_WITH_PREV:
      return "left mismatch coax with prev";
    case CTD_RCOAX_WITH_NEXT:
      return "right mismatch coax with next";
    case CTD_RCOAX_WITH_PREV:
      return "right mismatch coax with prev";
    case CTD_FCOAX_WITH_NEXT:
      return "flush coax with next";
    case CTD_FCOAX_WITH_PREV:
      return "flush coax with prev";
    default:
      verify_expr(false, "bug");
  }
}

std::vector<std::string> Structure::Description(int nesting) const {
  std::vector<std::string> desc;
  desc.push_back(sfmt("%d - %s", nesting, ShortDesc().c_str()));
  for (const auto& note : notes)
    desc.push_back(" | " + note);
  for (int i = 0; i < int(branches.size()); ++i)
    desc.push_back(" |-- " + BranchDesc(i));
  for (const auto& branch : branches) {
    auto branch_desc = branch->Description(nesting + 1);
    desc.insert(desc.end(), branch_desc.begin(), branch_desc.end());
  }
  return desc;
}

void Structure::AddNote(std::string note, ...) {
  va_list l;
  va_start(l, note);
  notes.push_back(vsfmt(note.c_str(), l));
  va_end(l);
}

std::string HairpinLoopStructure::ShortDesc() const {
  return sfmt("Hairpin(%d, %d) - %de:%de", st, en, GetTotalEnergy(), GetSelfEnergy());
}

std::string InternalLoopStructure::ShortDesc() const {
  return sfmt("InternalLoop(%d, %d, %d, %d) - %de:%de",
      ost, oen, ist, ien, GetTotalEnergy(), GetSelfEnergy());
}

std::string MultiLoopStructure::ShortDesc() const {
  return sfmt("MultiLoop(%d, %d) - %de:%de", st, en, GetTotalEnergy(), GetSelfEnergy());
}

std::string StackingStructure::ShortDesc() const {
  return sfmt("Stacking(%d, %d) - %de:%de", st, en, GetTotalEnergy(), GetSelfEnergy());
}
}
}
