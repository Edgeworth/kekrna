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
      // ..)...(..
      //  en   st
      const int lspace = en, rspace = N - st - 1;
      const base_t stb = gr[st],
          st1b = rspace ? gr[st + 1] : base_t(-1),
          st2b = rspace > 1 ? gr[st + 2] : base_t(-1),
          enb = gr[en],
          en1b = lspace ? gr[en - 1] : base_t(-1),
          en2b = lspace > 1 ? gr[en - 2] : base_t(-1);

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
        if (lspace) {
          p += gpt[0][en - 1][PT_U];  // Left empty, right filled.
          p += gpt[st + 1][N - 1][PT_U] * gpt[0][en - 1][PT_U];  // Both filled.
        }

        const auto base_branch_cost = Boltzmann(gpc.augubranch[stb][enb] + gem.multiloop_hack_a);
        // Enclosing cases:
        // |   >)   (<   |
        if (lspace && rspace) p += base_branch_cost * gpt[st + 1][en - 1][PT_U2];
        // |   >)   (3<  | 3'
        if (lspace && rspace > 1)
          p += base_branch_cost * gpt[st + 2][en - 1][PT_U2] *
              Boltzmann(gem.dangle3[stb][st1b][enb]);
        // |  >5)   (<   | 5'
        if (lspace > 1 && rspace)
          p += base_branch_cost * gpt[st + 1][en - 2][PT_U2] *
              Boltzmann(gem.dangle5[stb][en1b][enb]);
        // |  >m)   (m<  | Terminal mismatch
        if (lspace > 1 && rspace > 1)
          p += base_branch_cost * gpt[st + 2][en - 2][PT_U2] *
              Boltzmann(gem.terminal[stb][st1b][en1b][enb]);

        // For right cases: TODO update desc strings
        // |  ).)   (.<   >(  | Right outer coax - exterior loop
        // |  ).)   (<   >.(  | Right left coax - exterior loop
        // |  ))   (<   >(  | Right flush coax - exterior loop

        const auto outer_coax = lspace && rspace ? gem.MismatchCoaxial(stb, st1b, en1b, enb) : 0;
        // Left cases:
        // Exterior loops:
        // No space after piv; Can form loop up to piv
        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < N; ++piv) {
          base_t pl1b = gr[piv - 1], plb = gr[piv];
          // |<   >)   (.(   ).<   >| Left right coax - exterior loop
          // form piv loop
          const auto left = gpt[piv + 1][N - 1][PT_U] + 1.0;

          if (lspace) {
            // |<   >.)   (.(   )<   >| Left outer coax - exterior loop
            // lspace > 0 && form piv loop
            const auto exteriorr2 = left * (lspace > 1 ? gpt[0][en - 2][PT_U] + 1.0 : 1.0);
            p += gpt[st + 2][piv][PT_P] * exteriorr2 * Boltzmann(gem.AuGuPenalty(stb, enb) +
                gem.AuGuPenalty(st2b, plb) + outer_coax);
          }

          const auto exteriorr1 =  left * (lspace ? gpt[0][en - 1][PT_U] + 1.0 : 1.0);
          p += gpt[st + 2][piv - 1][PT_P] * exteriorr1 * Boltzmann(gem.AuGuPenalty(stb, enb) +
              gem.AuGuPenalty(st2b, pl1b) + gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
          // |<   >)   ((   )<   >| Left flush coax - exterior loop
          // form piv loop
          p += gpt[st + 1][piv][PT_P] * exteriorr1 * Boltzmann(gem.AuGuPenalty(stb, enb) +
              gem.AuGuPenalty(st1b, plb) + gem.stack[stb][st1b][plb][enb]);

        }

        // All enclosing cases require at least lspace > 0 to form an enclosing loop.
        // No space after piv; Can form loop up to piv; lspace > 0
        if (lspace) {
          // 1 space after piv; Can form loop up to piv; lspace > 0
          for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < N - 1; ++piv) {
            base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1];
            // |  >)   (.(   ).<  | Left right coax
            // lspace > 0 && piv + 1 < N && form piv loop
            p += base_branch_cost * gpt[st + 2][piv - 1][PT_P] * gpt[piv + 1][en - 1][PT_U] *
                Boltzmann(gpc.augubranch[st2b][pl1b] + gem.MismatchCoaxial(pl1b, plb, st1b, st2b));
            // |  >)   ((   )<  | Left flush coax
            // lspace > 0 && piv + 1 < N && form piv loop
            p += base_branch_cost * gpt[st + 1][piv][PT_P] * gpt[piv + 1][en - 1][PT_U] *
                Boltzmann(gpc.augubranch[st1b][plb] + gem.stack[stb][st1b][plb][enb]);
            // |  ))   (<   >(  | Right flush coax
            // lspace > 0 && piv + 1 < N && form piv loop
            p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 1][en - 1][PT_P] *
                Boltzmann(gpc.augubranch[prb][en1b] + gem.stack[stb][prb][en1b][enb]);
          }
        }

        // To form an enclosing loop, require lspace > 1.
        if (lspace > 1) {
          // 1 space after piv; Can form loop up to piv; lspace > 1
          for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < N - 1; ++piv) {
            base_t plb = gr[piv], prb = gr[piv + 1];

            if (piv < N - 2) {
              // |  ).)   (<   >.(  | Right left coax
              // lspace > 1 && piv + 2 < N && form piv loop
              base_t pr1b = gr[piv + 2];
              p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 2][en - 2][PT_P] *
                  Boltzmann(gpc.augubranch[pr1b][en2b] +
                      gem.MismatchCoaxial(en2b, en1b, prb, pr1b));
            }
            // |  >.)   (.(   )<  | Left outer coax
            // lspace > 1 && piv + 1 < N && form piv loop
            p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[piv + 1][en - 2][PT_U] *
                Boltzmann(gpc.augubranch[st2b][plb] + outer_coax);
            // |  ).)   (.<   >(  | Right outer coax
            // lspace > 1 && piv + 1 < N && form piv loop
            p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[piv + 1][en - 2][PT_P] *
                Boltzmann(gpc.augubranch[prb][en2b] + outer_coax);
          }
        }

        // Right cases:
        // Exterior loops:
        for (int piv = 0; piv < st; ++piv) {
          base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];

          if (lspace > 1 && rspace > 1) {
            //   >.)   (.(   )<   Left outer coax
            p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[piv + 1][en - 2][PT_U] *
                Boltzmann(gpc.augubranch[st2b][plb] + outer_coax);
            // Exterior loop case:
            p += gpt[st + 2][piv][PT_P] * (gpt[piv + 1][N - 1][PT_U] + 1.0) *
                (gpt[0][en - 2][PT_U] + 1.0) * Boltzmann(gem.AuGuPenalty(st2b, plb) + outer_coax);
            // >(   ).)   (.< Right outer coax
            p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[piv + 1][en - 2][PT_P] *
                Boltzmann(gpc.augubranch[prb][en2b] + outer_coax);
            // Exterior loop case:
            p += gpt[st + 2][piv][PT_U] * gpt[piv + 1][en - 2][PT_P] *
                Boltzmann(gpc.augubranch[prb][en2b] + outer_coax);
          }

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
}

}
}
}
