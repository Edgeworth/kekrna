#include <stack>
#include "fold.h"
#include "array.h"

namespace kekrna {
namespace fold {

enum {
  DP_P,  // For the paired array.
  DP_U,  // For the unpaired array.
  DP_U2, // Contains at least two branches.
  DP_U_WC,  // Unpaired but must start with a branch not involved in a CTD interaction that is not GU.
  DP_U_GU,  // Unpaired but must start with a branch not involved in a CTD interaction that is GU.
  DP_U_RCOAX,  // Unpaired but must start with a branch involved in a right coaxial stack - includes energy for it.
  DP_SIZE
};

enum {
  EXT,
  EXT_WC,
  EXT_GU,
  EXT_RCOAX,
  EXT_SIZE
};

#define UPDATE_CACHE(a, value) \
  do { \
    energy_t macro_upd_value_ = (value); \
    if (macro_upd_value_ < constants::CAP_E && macro_upd_value_ < arr[sz][st][a]) { \
      /*printf("Upd %d %d %d %d => %d " #value "\n", st, en, a, arr[sz][st][a], macro_upd_value_);*/ \
      arr[sz][st][a] = macro_upd_value_; \
    } \
  } while (0)

#define UPDATE_EXT(a, na, psz, pst, value) \
  do { \
    energy_t macro_upd_value_ = (value) + exterior[en + 1][na]; \
    if (macro_upd_value_ < constants::CAP_E && macro_upd_value_ < exterior[st][a]) { \
      /* printf("Ext st %d en %d a %d pst %d pen %d na %d %d => %d " #value "\n", st, en, a, pst, pst + psz - 1, na, exterior[st][a], macro_upd_value_); */ \
      exterior[st][a] = macro_upd_value_; \
      exterior_sts[st][a] = std::make_tuple(psz, pst, en + 1, na); \
    } \
  } while (0)

inline bool IsNotLonely(int st, int en) {
  return (en - st - 3 >= constants::HAIRPIN_MIN_SZ && CanPair(r[st + 1], r[en - 1])) ||
      (st > 0 && en < int(r.size() - 1) && CanPair(r[st - 1], r[en + 1]));
}

energy_t Fold(std::unique_ptr<structure::Structure>* s) {
  int N = int(r.size());
  assert(N > 0);
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, DP_SIZE> arr(r.size() + 1);
  assert(arr[N - 1][0][DP_U_WC] == constants::MAX_E);

  // TODO: check bounds
  // TODO: au/gu penalty?
  // sz includes st and en.
  static_assert(constants::HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int sz = constants::HAIRPIN_MIN_SZ + 2; sz <= N; ++sz) {
    int en = sz - 1;
    for (int st = 0; st < N - sz + 1; ++st) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

      // Update paired - only if can actually pair.
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en)) {
        for (int isz = 0; isz <= std::min(constants::TWOLOOP_MAX_SZ, sz - 4 - constants::HAIRPIN_MIN_SZ); ++isz) {
          for (int ist = st + 1; ist < st + isz + 2; ++ist) {
            int ien = en - (st + isz - ist) - 2;
            UPDATE_CACHE(DP_P, energy::TwoLoop(st, en, ist, ien) + arr[ien - ist + 1][ist][DP_P]);
          }
        }
        // Hairpin loops.
        UPDATE_CACHE(DP_P, energy::HairpinEnergy(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = energy::AuGuPenalty(st, en) + multiloop_hack_a + multiloop_hack_b;

        // No stacking case.
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 2][st + 1][DP_U2]);
        // (3<   ><   >) 3'
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 3][st + 2][DP_U2] + dangle3_e[stb][st1b][enb]);
        // (<   ><   >5) 5'
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 3][st + 1][DP_U2] + dangle5_e[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        UPDATE_CACHE(DP_P, base_branch_cost + arr[sz - 4][st + 2][DP_U2] + terminal_e[stb][st1b][en1b][enb]);

        for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz < sz - 3 - constants::HAIRPIN_MIN_SZ; ++lpivsz) {
          int rpivsz = sz - lpivsz - 2;
          // Paired coaxial stacking cases:
          base_t pl1b = r[st + lpivsz - 1], plb = r[st + lpivsz], prb = r[st + lpivsz + 1], pr1b = r[st + lpivsz + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // l_lABrCD =
          //   l => place a paired on the left
          //   skip A bases on left of left branch, B bases on right of left branch
          //   skip C bases on left of right branch, D bases on right of right branch.
          auto l_l00r00_cost = base_branch_cost + arr[lpivsz][st + 1][DP_P] + multiloop_hack_b +
              energy::AuGuPenalty(st + 1, st + lpivsz) + arr[rpivsz][st + 1 + lpivsz][DP_U];
          auto l_l10r01_cost = base_branch_cost + arr[lpivsz - 1][st + 2][DP_P] + multiloop_hack_b +
              energy::AuGuPenalty(st + 2, st + lpivsz) + arr[rpivsz - 1][st + 1 + lpivsz][DP_U];
          auto l_l11r00_cost = base_branch_cost + arr[lpivsz - 2][st + 2][DP_P] + multiloop_hack_b +
              energy::AuGuPenalty(st + 2, st + lpivsz - 1) + arr[rpivsz][st + 1 + lpivsz][DP_U];

          auto r_l00r00_cost = base_branch_cost + arr[lpivsz][st + 1][DP_U] + multiloop_hack_b +
              energy::AuGuPenalty(st + 1 + lpivsz, en - 1) + arr[rpivsz][st + 1 + lpivsz][DP_P];
          auto r_l10r01_cost = base_branch_cost + arr[lpivsz - 1][st + 2][DP_U] + multiloop_hack_b +
              energy::AuGuPenalty(st + 1 + lpivsz, en - 2) + arr[rpivsz - 1][st + 1 + lpivsz][DP_P];
          auto r_l00r11_cost = base_branch_cost + arr[lpivsz][st + 1][DP_U] + multiloop_hack_b +
              energy::AuGuPenalty(st + 2 + lpivsz, en - 2) + arr[rpivsz - 2][st + 2 + lpivsz][DP_P];

          // (.(   )   .) Left left coax - P
          auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
          UPDATE_CACHE(DP_P, l_l10r01_cost + outer_coax);
          // (.(   ).   ) Left right coax
          auto l_inner_coax = energy::MismatchMediatedCoaxialEnergy(pl1b, plb, st1b, st2b);
          UPDATE_CACHE(DP_P, l_l11r00_cost + l_inner_coax);
          // (   .(   ).) Right left coax
          auto r_inner_coax = energy::MismatchMediatedCoaxialEnergy(en2b, en1b, prb, pr1b);
          UPDATE_CACHE(DP_P, r_l00r11_cost + r_inner_coax);
          // (.   (   ).) Right right coax
          UPDATE_CACHE(DP_P, r_l10r01_cost + outer_coax);
          // ((   )   ) Left flush coax
          UPDATE_CACHE(DP_P, l_l00r00_cost + stacking_e[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          UPDATE_CACHE(DP_P, r_l00r00_cost + stacking_e[stb][prb][en1b][enb]);
        }
      }

      // Update unpaired.
      // Choose |st| to be unpaired.
      if (sz) {
        UPDATE_CACHE(DP_U, arr[sz - 1][st + 1][DP_U]);
        UPDATE_CACHE(DP_U2, arr[sz - 1][st + 1][DP_U2]);
      }
      // Pair here.
      UPDATE_CACHE(DP_U, arr[sz][st][DP_P] + multiloop_hack_b + energy::AuGuPenalty(st, en));
      for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz <= sz; ++lpivsz) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        int rpivsz = sz - lpivsz;
        auto pb = r[st + lpivsz - 1], pl1b = r[st + lpivsz - 2];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[lpivsz][st][DP_P] + energy::AuGuPenalty(st, st + lpivsz - 1) + multiloop_hack_b;
        auto base01 = arr[lpivsz - 1][st][DP_P] + energy::AuGuPenalty(st, st + lpivsz - 2) + multiloop_hack_b;
        auto base10 = arr[lpivsz - 1][st + 1][DP_P] + energy::AuGuPenalty(st + 1, st + lpivsz - 1) + multiloop_hack_b;
        auto base11 = arr[lpivsz - 2][st + 1][DP_P] + energy::AuGuPenalty(st + 1, st + lpivsz - 2) + multiloop_hack_b;
        // Min is for either placing another unpaired or leaving it as nothing.
        auto right_unpaired = std::min(arr[rpivsz][st + lpivsz][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        UPDATE_CACHE(DP_U2, base00 + arr[rpivsz][st + lpivsz][DP_U]);
        auto val = base00 + right_unpaired;
        UPDATE_CACHE(DP_U, val);
        if (IsGu(stb, pb))
          UPDATE_CACHE(DP_U_GU, val);
        else
          UPDATE_CACHE(DP_U_WC, val);

        // (   )3<   > 3' - U
        UPDATE_CACHE(DP_U, base01 + dangle3_e[pl1b][pb][stb] + right_unpaired);
        UPDATE_CACHE(DP_U2, base01 + dangle3_e[pl1b][pb][stb] + arr[rpivsz][st + lpivsz][DP_U]);
        // 5(   )<   > 5' - U
        UPDATE_CACHE(DP_U, base10 + dangle5_e[pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base10 + dangle5_e[pb][stb][st1b] + arr[rpivsz][st + lpivsz][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        UPDATE_CACHE(DP_U, base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired);
        UPDATE_CACHE(DP_U2, base11 + terminal_e[pl1b][pb][stb][st1b] + arr[rpivsz][st + lpivsz][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + energy::MismatchMediatedCoaxialEnergy(pl1b, pb, stb, st1b) +
            std::min(arr[rpivsz][st + lpivsz][DP_U_WC], arr[rpivsz][st + lpivsz][DP_U_GU]);
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);

        // (   ).<(   ). > Right coax forward and backward
        val = base01 + arr[rpivsz][st + lpivsz][DP_U_RCOAX];
        UPDATE_CACHE(DP_U, val);
        UPDATE_CACHE(DP_U2, val);
        if (st > 0)
          UPDATE_CACHE(DP_U_RCOAX, base01 + energy::MismatchMediatedCoaxialEnergy(
              pl1b, pb, r[st - 1], stb) + right_unpaired);

        // There has to be remaining bases to even have a chance at these cases.
        if (rpivsz > 0) {
          auto pr1b = r[st + lpivsz];
          // (   )<(   ) > Flush coax - U
          val = base00 + stacking_e[pb][pr1b][pr1b ^ 3][stb] + arr[rpivsz][st + lpivsz][DP_U_WC];
          UPDATE_CACHE(DP_U, val);
          UPDATE_CACHE(DP_U2, val);
          if (pr1b == G || pr1b == U) {
            val = base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[rpivsz][st + lpivsz][DP_U_GU];
            UPDATE_CACHE(DP_U, val);
            UPDATE_CACHE(DP_U2, val);
          }
        }
      }

      en++;
    }
  }

  // Exterior loop calculation. There can be no paired base on exterior[en].
  array2d_t<energy_t, EXT_SIZE> exterior(std::size_t(N + 1));
  exterior[N][EXT] = 0;
  // Holds: sz, st, nst, na - nst and na for next index into itself; sz, st for the paired.
  array2d_t<std::tuple<int, int, int, int>, EXT_SIZE> exterior_sts(std::size_t(N), 0xFF);
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    exterior[st][EXT] = exterior[st + 1][EXT];
    for (int sz = constants::HAIRPIN_MIN_SZ + 2; sz < N - st + 1; ++sz) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      int en = st + sz - 1, rem = N - st - sz;
      auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];

      auto base00 = arr[sz][st][DP_P] + energy::AuGuPenalty(st, en);
      auto base01 = arr[sz - 1][st][DP_P] + energy::AuGuPenalty(st, en - 1);
      auto base10 = arr[sz - 1][st + 1][DP_P] + energy::AuGuPenalty(st + 1, en);
      auto base11 = arr[sz - 2][st + 1][DP_P] + energy::AuGuPenalty(st + 1, en - 1);

      // (   )<   >
      UPDATE_EXT(EXT, EXT, sz, st, base00);
      if (IsGu(stb, enb))
        UPDATE_EXT(EXT_GU, EXT, sz, st, base00);
      else
        UPDATE_EXT(EXT_WC, EXT, sz, st, base00);

      // (   )3<   > 3'
      UPDATE_EXT(EXT, EXT, sz - 1, st, base01 + dangle3_e[en1b][enb][stb]);
      // 5(   )<   > 5'
      UPDATE_EXT(EXT, EXT, sz - 1, st + 1, base10 + dangle5_e[enb][stb][st1b]);
      // .(   ).<   > Terminal mismatch
      UPDATE_EXT(EXT, EXT, sz - 2, st + 1, base11 + terminal_e[en1b][enb][stb][st1b]);
      // .(   ).<(   ) > Left coax  x
      auto val = base11 + energy::MismatchMediatedCoaxialEnergy(en1b, enb, stb, st1b);
      UPDATE_EXT(EXT, EXT_GU, sz - 2, st + 1, val);
      UPDATE_EXT(EXT, EXT_WC, sz - 2, st + 1, val);

      // (   ).<(   ). > Right coax forward
      UPDATE_EXT(EXT, EXT_RCOAX, sz - 1, st, base01);
      // (   ).<( * ). > Right coax backward
      if (st > 0)
        UPDATE_EXT(EXT_RCOAX, EXT, sz - 1, st, base01 + energy::MismatchMediatedCoaxialEnergy(
            en1b, enb, r[st - 1], stb));

      if (rem > 0) {
        // (   )<(   ) > Flush coax
        auto enrb = r[en + 1];
        UPDATE_EXT(EXT, EXT_WC, sz, st, base00 + stacking_e[enb][enrb][enrb ^ 3][stb]);
        if (enrb == G || enrb == U)
          UPDATE_EXT(EXT, EXT_GU, sz, st, base00 + stacking_e[enb][enrb][enrb ^ 1][stb]);
      }
    }
  }

  // Backtrace.
  p = std::vector<int>(r.size(), -1);
  // sz, st, paired
  std::stack<std::tuple<int, int, int>> q;
  int ext_st = 0, ext_a = EXT;
  while (ext_st < N) {
    int psz, pst, nst, na;
    std::tie(psz, pst, nst, na) = exterior_sts[ext_st][ext_a];
    if (nst == -1) {
      ++ext_st;
      continue;
    }
    //printf("Exterior: %d %d %d %d\n", psz, pst, nst, na);
    q.emplace(psz, pst, DP_P);
    ext_st = nst;
    ext_a = na;
  }

  while (!q.empty()) {
    int sz, st, a;
    std::tie(sz, st, a) = q.top();
    int en = st + sz - 1;
    assert(sz >= constants::HAIRPIN_MIN_SZ + 2);
    auto stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];
    // printf("%d %d %d: %de\n", st, en, a, arr[sz][st][a]);
    q.pop();
    if (a == DP_P) {
      // It's paired, so add it to the folding.
      p[st] = en;
      p[en] = st;

      // Following largely matches the above DP so look up there for comments.
      for (int isz = 0; isz <= std::min(constants::TWOLOOP_MAX_SZ, sz - 4 - constants::HAIRPIN_MIN_SZ); ++isz) {
        for (int ist = st + 1; ist < st + isz + 2; ++ist) {
          int ien = en - (st + isz - ist) - 2;
          auto val = energy::TwoLoop(st, en, ist, ien) + arr[sz - isz - 2][ist][DP_P];
          if (val == arr[sz][st][DP_P]) {
            q.emplace(sz - isz - 2, ist, DP_P);
            goto loopend;
          }
        }
      }

      auto base_branch_cost = energy::AuGuPenalty(st, en) + multiloop_hack_a + multiloop_hack_b;

      // (<   ><    >)
      if (base_branch_cost + arr[sz - 2][st + 1][DP_U2] == arr[sz][st][DP_P]) {
        q.emplace(sz - 2, st + 1, DP_U2);
        goto loopend;
      }
      // (3<   ><   >) 3'
      if (base_branch_cost + arr[sz - 3][st + 2][DP_U2] + dangle3_e[stb][st1b][enb] == arr[sz][st][DP_P]) {
        q.emplace(sz - 3, st + 2, DP_U2);
        goto loopend;
      }
      // (<   ><   >5) 5'
      if (base_branch_cost + arr[sz - 3][st + 1][DP_U2] + dangle5_e[stb][en1b][enb] == arr[sz][st][DP_P]) {
        q.emplace(sz - 3, st + 1, DP_U2);
        goto loopend;
      }
      // (.<   ><   >.) Terminal mismatch
      if (base_branch_cost + arr[sz - 4][st + 2][DP_U2] + terminal_e[stb][st1b][en1b][enb] == arr[sz][st][DP_P]) {
        q.emplace(sz - 4, st + 2, DP_U2);
        goto loopend;
      }

      for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz < sz - 3 - constants::HAIRPIN_MIN_SZ; ++lpivsz) {
        auto rpivsz = sz - lpivsz - 2;
        base_t pl1b = r[st + lpivsz - 1], plb = r[st + lpivsz], prb = r[st + lpivsz + 1], pr1b = r[st + lpivsz + 2];

        auto l_l00r00_cost = base_branch_cost + arr[lpivsz][st + 1][DP_P] + multiloop_hack_b +
            energy::AuGuPenalty(st + 1, st + lpivsz) + arr[rpivsz][st + 1 + lpivsz][DP_U];
        auto l_l10r01_cost = base_branch_cost + arr[lpivsz - 1][st + 2][DP_P] + multiloop_hack_b +
            energy::AuGuPenalty(st + 2, st + lpivsz) + arr[rpivsz - 1][st + 1 + lpivsz][DP_U];
        auto l_l11r00_cost = base_branch_cost + arr[lpivsz - 2][st + 2][DP_P] + multiloop_hack_b +
            energy::AuGuPenalty(st + 2, st + lpivsz - 1) + arr[rpivsz][st + 1 + lpivsz][DP_U];

        auto r_l00r00_cost = base_branch_cost + arr[lpivsz][st + 1][DP_U] + multiloop_hack_b +
            energy::AuGuPenalty(st + 1 + lpivsz, en - 1) + arr[rpivsz][st + 1 + lpivsz][DP_P];
        auto r_l10r01_cost = base_branch_cost + arr[lpivsz - 1][st + 2][DP_U] + multiloop_hack_b +
            energy::AuGuPenalty(st + 1 + lpivsz, en - 2) + arr[rpivsz - 1][st + 1 + lpivsz][DP_P];
        auto r_l00r11_cost = base_branch_cost + arr[lpivsz][st + 1][DP_U] + multiloop_hack_b +
            energy::AuGuPenalty(st + 2 + lpivsz, en - 2) + arr[rpivsz - 2][st + 2 + lpivsz][DP_P];
        // (.(   )   .) Left left coax - P
        auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
        if (l_l10r01_cost + outer_coax == arr[sz][st][DP_P]) {
          q.emplace(lpivsz - 1, st + 2, DP_P);
          q.emplace(rpivsz - 1, st + 1 + lpivsz, DP_U);
          goto loopend;
        }
        // (.(   ).   ) Left right coax
        auto l_inner_coax = energy::MismatchMediatedCoaxialEnergy(pl1b, plb, st1b, st2b);
        if (l_l11r00_cost + l_inner_coax == arr[sz][st][DP_P]) {
          q.emplace(lpivsz - 2, st + 2, DP_P);
          q.emplace(rpivsz, st + 1 + lpivsz, DP_U);
          goto loopend;
        }
        // (   .(   ).) Right left coax
        auto r_inner_coax = energy::MismatchMediatedCoaxialEnergy(en2b, en1b, prb, pr1b);
        if (r_l00r11_cost + r_inner_coax == arr[sz][st][DP_P]) {
          q.emplace(lpivsz, st + 1, DP_U);
          q.emplace(rpivsz - 2, st + 2 + lpivsz, DP_P);
          goto loopend;
        }
        // (.   (   ).) Right right coax
        if (r_l10r01_cost + outer_coax == arr[sz][st][DP_P]) {
          q.emplace(lpivsz - 1, st + 2, DP_U);
          q.emplace(rpivsz - 1, st + 1 + lpivsz, DP_P);
          goto loopend;
        }
        // ((   )   ) Left flush coax
        if (l_l00r00_cost + stacking_e[stb][st1b][plb][enb] == arr[sz][st][DP_P]) {
          q.emplace(lpivsz, st + 1, DP_P);
          q.emplace(rpivsz, st + 1 + lpivsz, DP_U);
          goto loopend;
        }
        // (   (   )) Right flush coax
        if (r_l00r00_cost + stacking_e[stb][prb][en1b][enb] == arr[sz][st][DP_P]) {
          q.emplace(lpivsz, st + 1, DP_U);
          q.emplace(rpivsz, st + 1 + lpivsz, DP_P);
          goto loopend;
        }
      }
    } else {
      // Left unpaired. Either DP_U or DP_U2.
      if (sz && (a == DP_U || a == DP_U2) && arr[sz - 1][st + 1][a] == arr[sz][st][a]) {
        q.emplace(sz - 1, st + 1, a);
        goto loopend;
      }
      // Pair here. Only for DP_U.
      if (a == DP_U && arr[sz][st][DP_P] + multiloop_hack_b + energy::AuGuPenalty(st, en) == arr[sz][st][DP_U]) {
        q.emplace(sz, st, DP_P);
        goto loopend;
      }

      // Pair here.
      for (int lpivsz = constants::HAIRPIN_MIN_SZ + 2; lpivsz <= sz; ++lpivsz) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        int rpivsz = sz - lpivsz;
        auto pb = r[st + lpivsz - 1], pl1b = r[st + lpivsz - 2];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[lpivsz][st][DP_P] + energy::AuGuPenalty(st, st + lpivsz - 1) + multiloop_hack_b;
        auto base01 = arr[lpivsz - 1][st][DP_P] + energy::AuGuPenalty(st, st + lpivsz - 2) + multiloop_hack_b;
        auto base10 = arr[lpivsz - 1][st + 1][DP_P] + energy::AuGuPenalty(st + 1, st + lpivsz - 1) + multiloop_hack_b;
        auto base11 = arr[lpivsz - 2][st + 1][DP_P] + energy::AuGuPenalty(st + 1, st + lpivsz - 2) + multiloop_hack_b;

        // Min is for either placing another unpaired or leaving it as nothing.
        // If we're at U2, don't allow leaving as nothing.
        auto right_unpaired = arr[rpivsz][st + lpivsz][DP_U];
        if (a != DP_U2)
          right_unpaired = std::min(right_unpaired, 0);

        // Check a == U_RCOAX:
        // (   ).<( ** ). > Right coax backward
        if (a == DP_U_RCOAX) {
          if (st > 0 && base01 + energy::MismatchMediatedCoaxialEnergy(
              pl1b, pb, r[st - 1], stb) + right_unpaired == arr[sz][st][DP_U_RCOAX]) {
            q.emplace(lpivsz - 1, st, DP_P);
            if (right_unpaired)
              q.emplace(rpivsz, st + lpivsz, DP_U);
            goto loopend;

          }
          continue;
        }

        // (   )<   > - U, U2, U_WC?, U_GU?
        if (base00 + right_unpaired == arr[sz][st][a] &&
            (a != DP_U_WC || IsWatsonCrick(stb, pb)) &&
            (a != DP_U_GU || IsGu(stb, pb))) {
          q.emplace(lpivsz, st, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(rpivsz, st + lpivsz, DP_U);
          goto loopend;
        }

        // The rest of the cases are for U and U2.
        if (a != DP_U && a != DP_U2)
          continue;

        // (   )3<   > 3' - U, U2
        if (base01 + dangle3_e[pl1b][pb][stb] + right_unpaired == arr[sz][st][a]) {
          q.emplace(lpivsz - 1, st, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(rpivsz, st + lpivsz, DP_U);
          goto loopend;
        }
        // 5(   )<   > 5' - U, U2
        if (base10 + dangle5_e[pb][stb][st1b] + right_unpaired == arr[sz][st][a]) {
          q.emplace(lpivsz - 1, st + 1, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(rpivsz, st + lpivsz, DP_U);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch - U, U2
        if (base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired == arr[sz][st][a]) {
          q.emplace(lpivsz - 2, st + 1, DP_P);
          if (a == DP_U2 || right_unpaired)
            q.emplace(rpivsz, st + lpivsz, DP_U);
          goto loopend;
        }
        // .(   ).<(   ) > Left coax - U, U2
        auto val = base11 + energy::MismatchMediatedCoaxialEnergy(pl1b, pb, stb, st1b);
        if (val + arr[rpivsz][st + lpivsz][DP_U_WC] == arr[sz][st][a]) {
          q.emplace(lpivsz - 2, st + 1, DP_P);
          q.emplace(rpivsz, st + lpivsz, DP_U_WC);
          goto loopend;
        }
        if (val + arr[rpivsz][st + lpivsz][DP_U_GU] == arr[sz][st][a]) {
          q.emplace(lpivsz - 2, st + 1, DP_P);
          q.emplace(rpivsz, st + lpivsz, DP_U_GU);
          goto loopend;
        }

        // (   ).<(   ). > Right coax forward - U, U2
        if (base01 + arr[rpivsz][st + lpivsz][DP_U_RCOAX] == arr[sz][st][a]) {
          q.emplace(lpivsz - 1, st, DP_P);
          q.emplace(rpivsz, st + lpivsz, DP_U_RCOAX);
          goto loopend;
        }

        // There has to be remaining bases to even have a chance at these cases.
        if (rpivsz > 0) {
          auto pr1b = r[st + lpivsz];
          // (   )<(   ) > Flush coax - U, U2
          if (base00 + stacking_e[pb][pr1b][pr1b ^ 3][stb] + arr[rpivsz][st + lpivsz][DP_U_WC] == arr[sz][st][a]) {
            q.emplace(lpivsz, st, DP_P);
            q.emplace(rpivsz, st + lpivsz, DP_U_WC);
            goto loopend;
          }
          if ((pr1b == G || pr1b == U) &&
              base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[rpivsz][st + lpivsz][DP_U_GU] == arr[sz][st][a]) {
            q.emplace(lpivsz, st, DP_P);
            q.emplace(rpivsz, st + lpivsz, DP_U_GU);
            goto loopend;
          }
        }
      }
    }

    loopend:;
  }
  return exterior[0][EXT];
}

#undef UPDATE_CACHE
#undef UPDATE_EXT

namespace {

energy_t best;
std::vector<int> best_p;
std::vector<std::pair<int, int>> base_pairs;

void FoldBruteForceInternal(int idx) {
  if (idx == int(base_pairs.size())) {
    energy_t e = energy::ComputeEnergy(nullptr);
    if (e < best) {
      best = e;
      best_p = p;
    }
    return;
  }
  // Don't take this base pair.
  FoldBruteForceInternal(idx + 1);

  // Take this base pair.
  bool can_take = true;
  const auto& pair = base_pairs[idx];
  for (int i = pair.first; i <= pair.second; ++i) {
    if (p[i] != -1) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    p[pair.first] = pair.second;
    p[pair.second] = pair.first;
    FoldBruteForceInternal(idx + 1);
    p[pair.first] = -1;
    p[pair.second] = -1;
  }
}

}

energy_t FoldBruteForce(std::unique_ptr<structure::Structure>* s) {
  p = std::vector<int>(r.size(), -1);
  best = constants::MAX_E;
  base_pairs.clear();
  for (int st = 0; st < int(r.size()); ++st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < int(r.size()); ++en) {
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en))
        base_pairs.emplace_back(st, en);
    }
  }
  FoldBruteForceInternal(0);
  p = best_p;
  return energy::ComputeEnergy(s);
}

}
}
