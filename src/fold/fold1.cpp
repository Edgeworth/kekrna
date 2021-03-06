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
#include "fold/fold.h"

namespace kekrna {
namespace fold {
namespace internal {

using namespace energy;

void ComputeTables1() {
  const int N = int(gr.size());
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const base_t stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
          en1b = gr[en - 1], en2b = gr[en - 2];

      // Update paired - only if can actually pair.
      if (ViableFoldingPair(st, en)) {
        energy_t p_min = MAX_E;
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
            if (gdp[ist][ien][DP_P] < CAP_E)
              p_min = std::min(p_min, FastTwoLoop(st, en, ist, ien) + gdp[ist][ien][DP_P]);
          }
        }
        // Hairpin loops.
        p_min = std::min(p_min, gem.Hairpin(gr, st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost = gpc.augubranch[stb][enb] + gem.multiloop_hack_a;

        // (<   ><   >)
        p_min = std::min(p_min, base_branch_cost + gdp[st + 1][en - 1][DP_U2]);
        // (3<   ><   >) 3'
        p_min = std::min(
            p_min, base_branch_cost + gdp[st + 2][en - 1][DP_U2] + gem.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p_min = std::min(
            p_min, base_branch_cost + gdp[st + 1][en - 2][DP_U2] + gem.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p_min = std::min(p_min,
            base_branch_cost + gdp[st + 2][en - 2][DP_U2] + gem.terminal[stb][st1b][en1b][enb]);

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
          p_min = std::min(p_min, base_branch_cost + gdp[st + 2][piv][DP_P] +
              gpc.augubranch[st2b][plb] + gdp[piv + 1][en - 2][DP_U] + outer_coax);
          // (.   (   ).) Right outer coax
          p_min = std::min(p_min, base_branch_cost + gdp[st + 2][piv][DP_U] +
              gpc.augubranch[prb][en2b] + gdp[piv + 1][en - 2][DP_P] + outer_coax);

          // (.(   ).   ) Left right coax
          p_min = std::min(p_min, base_branch_cost + gdp[st + 2][piv - 1][DP_P] +
              gpc.augubranch[st2b][pl1b] + gdp[piv + 1][en - 1][DP_U] +
              gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          p_min = std::min(p_min, base_branch_cost + gdp[st + 1][piv][DP_U] +
              gpc.augubranch[pr1b][en2b] + gdp[piv + 2][en - 2][DP_P] +
              gem.MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p_min = std::min(p_min, base_branch_cost + gdp[st + 1][piv][DP_P] +
              gpc.augubranch[st1b][plb] + gdp[piv + 1][en - 1][DP_U] +
              gem.stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p_min = std::min(p_min, base_branch_cost + gdp[st + 1][piv][DP_U] +
              gpc.augubranch[prb][en1b] + gdp[piv + 1][en - 1][DP_P] +
              gem.stack[stb][prb][en1b][enb]);
        }

        gdp[st][en][DP_P] = p_min;
      }
      energy_t u_min = MAX_E, u2_min = MAX_E, rcoax_min = MAX_E, wc_min = MAX_E, gu_min = MAX_E;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u_min = std::min(u_min, gdp[st + 1][en][DP_U]);
        u2_min = std::min(u2_min, gdp[st + 1][en][DP_U2]);
      }
      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = gr[piv], pl1b = gr[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = gdp[st][piv][DP_P] + gpc.augubranch[stb][pb];
        const auto base01 = gdp[st][piv - 1][DP_P] + gpc.augubranch[stb][pl1b];
        const auto base10 = gdp[st + 1][piv][DP_P] + gpc.augubranch[st1b][pb];
        const auto base11 = gdp[st + 1][piv - 1][DP_P] + gpc.augubranch[st1b][pl1b];
        // Min is for either placing another unpaired or leaving it as nothing.
        const auto right_unpaired = std::min(gdp[piv + 1][en][DP_U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        u2_min = std::min(u2_min, base00 + gdp[piv + 1][en][DP_U]);
        auto val = base00 + right_unpaired;
        u_min = std::min(u_min, val);
        if (IsGu(stb, pb))
          gu_min = std::min(gu_min, val);
        else
          wc_min = std::min(wc_min, val);

        // (   )3<   > 3' - U
        u_min = std::min(u_min, base01 + gem.dangle3[pl1b][pb][stb] + right_unpaired);
        u2_min = std::min(u2_min, base01 + gem.dangle3[pl1b][pb][stb] + gdp[piv + 1][en][DP_U]);
        // 5(   )<   > 5' - U
        u_min = std::min(u_min, base10 + gem.dangle5[pb][stb][st1b] + right_unpaired);
        u2_min = std::min(u2_min, base10 + gem.dangle5[pb][stb][st1b] + gdp[piv + 1][en][DP_U]);
        // .(   ).<   > Terminal mismatch - U
        u_min = std::min(u_min, base11 + gem.terminal[pl1b][pb][stb][st1b] + right_unpaired);
        u2_min =
            std::min(u2_min, base11 + gem.terminal[pl1b][pb][stb][st1b] + gdp[piv + 1][en][DP_U]);
        // .(   ).<(   ) > Left coax - U
        val = base11 + gem.MismatchCoaxial(pl1b, pb, stb, st1b) +
            std::min(gdp[piv + 1][en][DP_U_WC], gdp[piv + 1][en][DP_U_GU]);
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);

        // (   ).<(   ). > Right coax forward and backward
        val = base01 + gdp[piv + 1][en][DP_U_RCOAX];
        u_min = std::min(u_min, val);
        u2_min = std::min(u2_min, val);
        if (st > 0)
          rcoax_min = std::min(
              rcoax_min, base01 + gem.MismatchCoaxial(pl1b, pb, gr[st - 1], stb) + right_unpaired);

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          const auto pr1b = gr[piv + 1];
          // (   )<(   ) > Flush coax - U
          val = base00 + gem.stack[pb][pr1b][pr1b ^ 3][stb] + gdp[piv + 1][en][DP_U_WC];
          u_min = std::min(u_min, val);
          u2_min = std::min(u2_min, val);
          if (pr1b == G || pr1b == U) {
            val = base00 + gem.stack[pb][pr1b][pr1b ^ 1][stb] + gdp[piv + 1][en][DP_U_GU];
            u_min = std::min(u_min, val);
            u2_min = std::min(u2_min, val);
          }
        }
      }

      gdp[st][en][DP_U] = u_min;
      gdp[st][en][DP_U2] = u2_min;
      gdp[st][en][DP_U_WC] = wc_min;
      gdp[st][en][DP_U_GU] = gu_min;
      gdp[st][en][DP_U_RCOAX] = rcoax_min;
    }
  }
}
}
}
}
