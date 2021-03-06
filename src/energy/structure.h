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
#ifndef KEKRNA_STRUCTURE_H
#define KEKRNA_STRUCTURE_H

#include <memory>
#include "common.h"
#include "energy/energy.h"
#include "energy/energy_internal.h"

namespace kekrna {
namespace energy {

const char* CtdToName(Ctd ctd);

class Structure {
public:
  Structure() = default;
  Structure(const Structure&) = delete;
  Structure& operator=(const Structure&) = delete;

  virtual ~Structure() = default;

  // This is not a reference because of varargs.
  void AddNote(std::string note, ...);
  std::vector<std::string> Description(int nesting = 0) const;
  virtual std::string ShortDesc() const = 0;

  virtual void AddBranch(std::unique_ptr<Structure> b) { branches.push_back(std::move(b)); }
  virtual std::string BranchDesc(int idx) const { return branches[idx]->ShortDesc(); }

  void SetSelfEnergy(energy_t e) { self_energy = e; }
  void SetTotalEnergy(energy_t e) { total_energy = e; }
  energy_t GetSelfEnergy() const { return self_energy; }
  energy_t GetTotalEnergy() const { return total_energy; }

protected:
  std::vector<std::unique_ptr<Structure>> branches;

private:
  energy_t self_energy;
  energy_t total_energy;
  std::vector<std::string> notes;
};

class HairpinLoopStructure : public Structure {
public:
  HairpinLoopStructure(int st_, int en_) : st(st_), en(en_) {}

  void AddBranch(std::unique_ptr<Structure>) { assert(false); }
  std::string ShortDesc() const;

private:
  int st, en;
};

class InternalLoopStructure : public Structure {
public:
  InternalLoopStructure(int ost_, int oen_, int ist_, int ien_)
      : ost(ost_), oen(oen_), ist(ist_), ien(ien_) {}

  void AddBranch(std::unique_ptr<Structure> b) {
    assert(branches.empty());
    Structure::AddBranch(std::move(b));
  }
  std::string ShortDesc() const;

private:
  int ost, oen, ist, ien;
};

class MultiLoopStructure : public Structure {
public:
  MultiLoopStructure(int st_, int en_) : st(st_), en(en_) {}

  void AddCtd(Ctd ctd, energy_t ctd_energy) { branch_ctds.emplace_back(ctd, ctd_energy); }
  std::string BranchDesc(int idx) const {
    return sfmt("%s - %de - %s", branches[idx]->ShortDesc().c_str(), branch_ctds[idx].second,
        CtdToName(branch_ctds[idx].first));
  }
  std::string ShortDesc() const;

private:
  int st, en;
  internal::branch_ctd_t branch_ctds;
};

class StackingStructure : public Structure {
public:
  StackingStructure(int st_, int en_) : st(st_), en(en_) {}

  void AddBranch(std::unique_ptr<Structure> b) {
    assert(branches.empty());
    Structure::AddBranch(std::move(b));
  }
  std::string ShortDesc() const;

private:
  int st, en;
};
}
}

#endif  // KEKRNA_PARSING_H
