#ifndef KEKRNA_FOLD_INTERNAL_H
#define KEKRNA_FOLD_INTERNAL_H

#include "common.h"
#include "base.h"
#include "fold/globals.h"

namespace kekrna {
namespace fold {
namespace internal {

inline bool ViableFoldingPair(int st, int en) {
  return CanPair(gr[st], gr[en]) &&
      ((en - st - 3 >= HAIRPIN_MIN_SZ && CanPair(gr[st + 1], gr[en - 1])) ||
          (st > 0 && en < int(gr.size() - 1) && CanPair(gr[st - 1], gr[en + 1])));
}

struct cand_t {
  energy_t energy;
  int idx;
};

struct index_t {
  int st, en, a;

  index_t() : st(-1), en(-1), a(-1) {}
  index_t(int st_, int en_, int a_) : st(st_), en(en_), a(a_) {}

  bool operator==(const index_t& o) const {
    return st == o.st && en == o.en && a == o.a;
  }

  bool operator!=(const index_t& o) const {
    return !(*this == o);
  }

  bool operator<(const index_t& o) const {
    if (st != o.st) return st < o.st;
    if (en != o.en) return en < o.en;
    return a < o.a;
  }
};

struct ctd_idx_t {
  ctd_idx_t() : idx(-1), ctd(CTD_NA) {}
  ctd_idx_t(int idx_, Ctd ctd_) : idx(idx_), ctd(ctd_) {}

  int idx;
  Ctd ctd;

  bool operator<(const ctd_idx_t& o) const {
    if (idx != o.idx) return idx < idx;
    return ctd < o.ctd;
  }

  bool operator==(const ctd_idx_t& o) const {
    return idx == o.idx && ctd == o.ctd;
  }

  bool operator!=(const ctd_idx_t& o) const {
    return !(*this == o);
  }
};

void ComputeTables0();
void ComputeTables1();
void ComputeTables2();
void ComputeTables3();
void ComputeExterior();
void Traceback();

}
}
}

namespace std {
template<>
struct hash<kekrna::fold::internal::index_t> {
  size_t operator()(const kekrna::fold::internal::index_t& o) const {
    size_t v = 0;
    v ^= hash<int>()(o.st) + 0x9e3779b9 + (v << 6) + (v >> 2);
    v ^= hash<int>()(o.en) + 0x9e3779b9 + (v << 6) + (v >> 2);
    v ^= hash<int>()(o.a) + 0x9e3779b9 + (v << 6) + (v >> 2);
    return v;
  }
};
}

#endif //KEKRNA_FOLD_INTERNAL_H
