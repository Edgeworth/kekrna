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
#include "partition/partition.h"
#include "globals.h"
#include "partition/partition_globals.h"
#include "energy/energy_globals.h"

namespace kekrna {
namespace partition {
namespace internal {

using energy::gem;
using energy::gpc;
using energy::Boltzmann;

void Partition0() {
  // TODO do enclosing stuff
  const int N = int(gr.size());
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const base_t stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
          en1b = gr[en - 1], en2b = gr[en - 2];

      //if (CanPair(stb, enb)) {  // TODO lonely pairs?
      if (energy::ViableFoldingPair(st, en)) {
        penergy_t p = 0.0;
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist)  // TODO not accurate
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien)
            p += Boltzmann(energy::FastTwoLoop(st, en, ist, ien)) * gpt[ist][ien][PT_P];
        // Hairpin loops.
        p += Boltzmann(gem.Hairpin(gr, st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const auto base_branch_cost = Boltzmann(gpc.augubranch[stb][enb] + gem.multiloop_hack_a);

        // (<   ><   >)
        p += base_branch_cost * gpt[st + 1][en - 1][PT_U2];
        // (3<   ><   >) 3'
        p += base_branch_cost * gpt[st + 2][en - 1][PT_U2] * Boltzmann(gem.dangle3[stb][st1b][enb]);
        // (<   ><   >5) 5'
        p += base_branch_cost * gpt[st + 1][en - 2][PT_U2] * Boltzmann(gem.dangle5[stb][en1b][enb]);
        // (.<   ><   >.) Terminal mismatch
        p += base_branch_cost * gpt[st + 2][en - 2][PT_U2] *
            Boltzmann(gem.terminal[stb][st1b][en1b][enb]);

        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];
          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // (.(   )   .) Left outer coax - P
          const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);
          p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[piv + 1][en - 2][PT_U] *
              Boltzmann(gpc.augubranch[st2b][plb] + outer_coax);
          // (.   (   ).) Right outer coax
          p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[piv + 1][en - 2][PT_P] *
              Boltzmann(gpc.augubranch[prb][en2b] + outer_coax);

          // (.(   ).   ) Left right coax
          p += base_branch_cost * gpt[st + 2][piv - 1][PT_P] * gpt[piv + 1][en - 1][PT_U] *
              Boltzmann(gpc.augubranch[st2b][pl1b] + gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
          // (   .(   ).) Right left coax
          p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 2][en - 2][PT_P] *
              Boltzmann(gpc.augubranch[pr1b][en2b] + gem.MismatchCoaxial(en2b, en1b, prb, pr1b));

          // ((   )   ) Left flush coax
          p += base_branch_cost * gpt[st + 1][piv][PT_P] * gpt[piv + 1][en - 1][PT_U] *
              Boltzmann(gpc.augubranch[st1b][plb] + gem.stack[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 1][en - 1][PT_P] *
              Boltzmann(gpc.augubranch[prb][en1b] + gem.stack[stb][prb][en1b][enb]);
        }

        gpt[st][en][PT_P] = p;
      }
      penergy_t u = 0.0, u2 = 0.0, rcoax = 0.0, wc = 0.0, gu = 0.0;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u += gpt[st + 1][en][PT_U];
        u2 += gpt[st + 1][en][PT_U2];
      }

      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = gr[piv], pl1b = gr[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = gpt[st][piv][PT_P] * Boltzmann(gpc.augubranch[stb][pb]);
        const auto base01 = gpt[st][piv - 1][PT_P] * Boltzmann(gpc.augubranch[stb][pl1b]);
        const auto base10 = gpt[st + 1][piv][PT_P] * Boltzmann(gpc.augubranch[st1b][pb]);
        const auto base11 = gpt[st + 1][piv - 1][PT_P] * Boltzmann(gpc.augubranch[st1b][pl1b]);

        // (   )<   > - U, U_WC?, U_GU?
        u2 += base00 * gpt[piv +
            1][en][PT_U];  // TODO accesses like this when add stuff in other part of table
        auto val = base00 + base00 * gpt[piv + 1][en][PT_U];
        u += val;
        if (IsGu(stb, pb)) gu += val;
        else wc += val;

        // (   )3<   > 3' - U
        val = base01 * Boltzmann(gem.dangle3[pl1b][pb][stb]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // 5(   )<   > 5' - U
        val = base10 * Boltzmann(gem.dangle5[pb][stb][st1b]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<   > Terminal mismatch - U
        val = base11 * Boltzmann(gem.terminal[pl1b][pb][stb][st1b]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<(   ) > Left coax - U
        val = base11 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, stb, st1b));
        val = val * (gpt[piv + 1][en][PT_U_WC] + gpt[piv + 1][en][PT_U_GU]);
        u += val;
        u2 += val;

        // (   ).<(   ). > Right coax forward and backward
        val = base01 * gpt[piv + 1][en][PT_U_RCOAX];
        u += val;
        u2 += val;
        if (st > 0) {
          val = base01 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, gr[st - 1], stb));
          rcoax += val;
          rcoax += val * gpt[piv + 1][en][PT_U];
        }

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          const auto pr1b = gr[piv + 1];
          // (   )<(   ) > Flush coax - U
          val = base00 * Boltzmann(gem.stack[pb][pr1b][pr1b ^ 3][stb]) * gpt[piv + 1][en][PT_U_WC];
          u += val;
          u2 += val;
          if (pr1b == G || pr1b == U) {
            val =
                base00 * Boltzmann(gem.stack[pb][pr1b][pr1b ^ 1][stb]) * gpt[piv + 1][en][PT_U_GU];
            u += val;
            u2 += val;
          }
        }
      }
      gpt[st][en][PT_U] = u;
      gpt[st][en][PT_U2] = u2;
      gpt[st][en][PT_U_WC] = wc;
      gpt[st][en][PT_U_GU] = gu;
      gpt[st][en][PT_U_RCOAX] = rcoax;
    }
  }

  // Fill the left triangle.
  // The meaning of the tables changes here:
  // U, U2: any table index with en < st must have a loop enclosing (en, st)
  for (int st = 0; st < N; ++st) {
    for (int en = 0; en < st; ++en) {
      //        ..)...(..
      // rspace  en   st  lspace
      const int lspace = N - st - 1, rspace = en;
      const base_t stb = gr[st],
          st1b = lspace ? gr[st + 1] : base_t(-1),
          st2b = lspace > 1 ? gr[st + 2] : base_t(-1),
          enb = gr[en],
          en1b = rspace ? gr[en - 1] : base_t(-1),
          en2b = rspace > 1 ? gr[en - 2] : base_t(-1);

      //if (CanPair(enb, stb)) {  // TODO lonely pairs?
      if (energy::ViableFoldingPair(en, st)) {
        penergy_t p = 0.0;
        const int ost_max = std::min(st + TWOLOOP_MAX_SZ + 2, N);
        for (int ost = st + 1; ost < ost_max; ++ost) {
          const int oen_min = std::min(en - TWOLOOP_MAX_SZ - 2 + (ost - st - 1), 0);
          for (int oen = en - 1; oen > oen_min; --oen)
            p += Boltzmann(energy::FastTwoLoop(oen, ost, en, st)) * gpt[ost][oen][PT_P];
        }

        // Try being an exterior loop - coax cases handled in the loop below.
        // TODO try to merge all the conditionals here
        p += Boltzmann(gem.AuGuPenalty(enb, stb));
        p += gpt[st + 1][N - 1][PT_U];  // Left filled, right empty.
        if (rspace) {
          p += gpt[0][en - 1][PT_U];  // Left empty, right filled.
          p += gpt[st + 1][N - 1][PT_U] * gpt[0][en - 1][PT_U];  // Both filled.
        }

        const auto base_branch_cost = Boltzmann(gpc.augubranch[stb][enb] + gem.multiloop_hack_a);
        // Enclosing cases:
        // |   >)   (<   |
        if (lspace && rspace) p += base_branch_cost * gpt[st + 1][en - 1][PT_U2];
        // |   >)   (3<  | 3'
        if (lspace > 1 && rspace)
          p += base_branch_cost * gpt[st + 2][en - 1][PT_U2] *
              Boltzmann(gem.dangle3[stb][st1b][enb]);
        // |  >5)   (<   | 5'
        if (lspace && rspace > 1)
          p += base_branch_cost * gpt[st + 1][en - 2][PT_U2] *
              Boltzmann(gem.dangle5[stb][en1b][enb]);
        // |  >m)   (m<  | Terminal mismatch
        if (lspace > 1 && rspace > 1)
          p += base_branch_cost * gpt[st + 2][en - 2][PT_U2] *
              Boltzmann(gem.terminal[stb][st1b][en1b][enb]);

        // TODO merge if statements?
        const int limit = en - HAIRPIN_MIN_SZ - 2 + (en < st ? N : 0);
        for (int tpiv = st + HAIRPIN_MIN_SZ + 2; tpiv < limit; ++tpiv) {
          const int pl = FastMod(tpiv - 1, N), piv = FastMod(tpiv, N),
              pr = FastMod(tpiv + 1, N), pr1 = FastMod(tpiv + 2, N);
          // Left block is: [st, piv], Right block is: [piv + 1, en].
          base_t pl1b = gr[pl], plb = gr[tpiv], prb = gr[pr], pr1b = gr[pr1];

          // When neither the left block nor the right block straddles the border
          // we don't get enclosing loops sometimes. So check that with |not_straddling|.
          const bool straddling = tpiv != (N - 1);
          const bool left_dot_straddling = straddling && tpiv != 0;  // Can't split a dot.
          const bool right_dot_straddling = straddling && tpiv != N - 2;  // Can't split a dot.
          const bool left_not_straddling = tpiv < N;  // Implies lspace > 1
          const bool right_not_straddling = tpiv >= N - 1;  // Implies rspace > 1

          if (lspace > 1 && rspace > 1 && straddling) {
            const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);  // TODO recomputation here
            // |  >.)   (.(   )<  | Enclosing loop - Left outer coax
            // lspace > 1 && rspace > 1 && enclosed
            p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[pr][en - 2][PT_U] *
                Boltzmann(gpc.augubranch[st2b][plb] + outer_coax);
            // |  >(   ).)   (.<  | Enclosing loop - Right outer coax
            // lspace > 1 && rspace > 1 && enclosed
            p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[pr][en - 2][PT_P] *
                Boltzmann(gpc.augubranch[prb][en2b] + outer_coax);
          }

          if (lspace && right_not_straddling) {
            const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);  // TODO recomputation here
            const auto left_exterior = lspace > 1 ? gpt[st + 2][N - 1][PT_U] + 1.0 : 1.0;
            const auto right_exterior = gpt[0][piv][PT_U] + 1.0;
            // |<   >(   ).)   (.<   >| Exterior loop - Right outer coax
            // lspace > 0 && rspace > 1 && not enclosed
            p += gpt[pr][en - 2][PT_P] * left_exterior * right_exterior * Boltzmann(
                gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(prb, en2b) + outer_coax);
          }

          if (rspace && left_not_straddling) {
            const auto outer_coax = gem.MismatchCoaxial(stb, st1b, en1b, enb);  // TODO recomputation here
            const auto left_exterior = gpt[piv + 1][N - 1][PT_U] + 1.0;
            const auto right_exterior = rspace > 1 ? gpt[0][en - 2][PT_U] + 1.0 : 1.0;
            // |<   >.)   (.(   )<   >| Exterior loop - Left outer coax
            // lspace > 0 && rspace > 0 && not enclosed
            p += gpt[st + 2][piv][PT_P] * left_exterior * right_exterior * Boltzmann(
                gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(st2b, plb) + outer_coax);
          }

          if (lspace > 1 && rspace && left_dot_straddling) {
            // |  >)   (.(   ).<  | Enclosing loop - Left right coax
            // lspace > 1 && rspace > 0 && enclosed && no dot split
            p += base_branch_cost * gpt[st + 2][pl][PT_P] * gpt[pr][en - 1][PT_U] *
                Boltzmann(gpc.augubranch[st2b][pl1b] + gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
          }

          if (left_not_straddling) {
            const auto left_exterior = gpt[piv + 1][N - 1][PT_U] + 1.0;
            const auto right_exterior = rspace ? gpt[0][en - 1][PT_U] + 1.0 : 1.0;
            // |<   >)   (.(   ).<   >| Exterior loop - Left right coax
            // lspace > 1 && not enclosed
            p += gpt[st + 2][pl][PT_P] * left_exterior * right_exterior * Boltzmann(
                gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(st2b, pl1b) +
                    gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
            // |<   >)   ((   )<   >| Exterior loop - Left flush coax
            // lspace > 0 && not enclosed
            p += gpt[st + 1][piv][PT_P] * left_exterior * right_exterior * Boltzmann(
                gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(st1b, plb) +
                    gem.stack[stb][st1b][plb][enb]);
          }

          if (lspace && rspace > 1 && right_dot_straddling) {
            // |  >.(   ).)   (<  | Enclosing loop - Right left coax
            // lspace > 0 && rspace > 1 && enclosed && no dot split
            p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[pr1][en - 2][PT_P] *
                Boltzmann(gpc.augubranch[pr1b][en2b] + gem.MismatchCoaxial(en2b, en1b, prb, pr1b));
          }

          if (right_not_straddling) {
            const auto left_exterior = gpt[st + 1][N - 1][PT_U] + 1.0;
            const auto right_exterior = gpt[0][piv][PT_U] + 1.0;
            // |<   >.(   ).)   (<   >| Exterior loop - Right left coax
            // rspace > 1 && not enclosed
            p += gpt[pr1][en - 2][PT_P] * left_exterior * right_exterior * Boltzmann(
                gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(pr1b, en2b) +
                    gem.MismatchCoaxial(en2b, en1b, plb, pr1b));
            // |<   >(   ))   (<   >| Exterior loop - Right flush coax
            // rspace > 0 && not enclosed
            p += gpt[pr][en - 1][PT_P] * left_exterior * right_exterior * Boltzmann(
                gem.AuGuPenalty(stb, enb) + gem.AuGuPenalty(prb, en1b) +
                    gem.stack[stb][prb][en1b][enb]);
          }

          if (lspace && rspace && straddling) {
            // |  >)   ((   )<  | Enclosing loop - Left flush coax
            // lspace > 0 && rspace > 0 && enclosed
            p += base_branch_cost * gpt[st + 1][piv][PT_P] * gpt[pr][en - 1][PT_U] *
                Boltzmann(gpc.augubranch[st1b][plb] + gem.stack[stb][st1b][plb][enb]);
            // |  >(   ))   (<  | Enclosing loop - Right flush coax
            // lspace > 0 && rspace > 0 && enclosed
            p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[pr][en - 1][PT_P] *
                Boltzmann(gpc.augubranch[prb][en1b] + gem.stack[stb][prb][en1b][enb]);
          }
        }

        gpt[st][en][PT_P] = p;
      }
      penergy_t u = 0.0, u2 = 0.0, rcoax = 0.0, wc = 0.0, gu = 0.0;
      // Update unpaired.
      // Choose |st| to be unpaired.
      if (st + 1 < en) {
        u += gpt[st + 1][en][PT_U];
        u2 += gpt[st + 1][en][PT_U2];
      }

      for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        const auto pb = gr[piv], pl1b = gr[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        const auto base00 = gpt[st][piv][PT_P] * Boltzmann(gpc.augubranch[stb][pb]);
        const auto base01 = gpt[st][piv - 1][PT_P] * Boltzmann(gpc.augubranch[stb][pl1b]);
        const auto base10 = gpt[st + 1][piv][PT_P] * Boltzmann(gpc.augubranch[st1b][pb]);
        const auto base11 = gpt[st + 1][piv - 1][PT_P] * Boltzmann(gpc.augubranch[st1b][pl1b]);

        // (   )<   > - U, U_WC?, U_GU?
        u2 += base00 * gpt[piv +
            1][en][PT_U];  // TODO accesses like this when add stuff in other part of table
        auto val = base00 + base00 * gpt[piv + 1][en][PT_U];
        u += val;
        if (IsGu(stb, pb)) gu += val;
        else wc += val;

        // (   )3<   > 3' - U
        val = base01 * Boltzmann(gem.dangle3[pl1b][pb][stb]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // 5(   )<   > 5' - U
        val = base10 * Boltzmann(gem.dangle5[pb][stb][st1b]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<   > Terminal mismatch - U
        val = base11 * Boltzmann(gem.terminal[pl1b][pb][stb][st1b]);
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<(   ) > Left coax - U
        val = base11 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, stb, st1b));
        val = val * (gpt[piv + 1][en][PT_U_WC] + gpt[piv + 1][en][PT_U_GU]);
        u += val;
        u2 += val;

        // (   ).<(   ). > Right coax forward and backward
        val = base01 * gpt[piv + 1][en][PT_U_RCOAX];
        u += val;
        u2 += val;
        if (st > 0) {
          val = base01 * Boltzmann(gem.MismatchCoaxial(pl1b, pb, gr[st - 1], stb));
          rcoax += val;
          rcoax += val * gpt[piv + 1][en][PT_U];
        }

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          const auto pr1b = gr[piv + 1];
          // (   )<(   ) > Flush coax - U
          val = base00 * Boltzmann(gem.stack[pb][pr1b][pr1b ^ 3][stb]) * gpt[piv + 1][en][PT_U_WC];
          u += val;
          u2 += val;
          if (pr1b == G || pr1b == U) {
            val =
                base00 * Boltzmann(gem.stack[pb][pr1b][pr1b ^ 1][stb]) * gpt[piv + 1][en][PT_U_GU];
            u += val;
            u2 += val;
          }
        }
      }
      gpt[st][en][PT_U] = u;
      gpt[st][en][PT_U2] = u2;
      gpt[st][en][PT_U_WC] = wc;
      gpt[st][en][PT_U_GU] = gu;
      gpt[st][en][PT_U_RCOAX] = rcoax;
    }
  }
}

}
}
}
