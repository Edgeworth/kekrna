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
#include <iostream>
#include <energy/energy_globals.h>
#include "partition/partition.h"
#include "globals.h"
#include "partition/partition_globals.h"
#include "energy/energy_globals.h"

namespace kekrna {
namespace partition {
namespace internal {

void Partition1() {
  const int N = int(gr.size());
  static_assert(
      HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int st = N - 1; st >= 0; --st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      const base_t stb = gr[st], st1b = gr[st + 1], st2b = gr[st + 2], enb = gr[en],
          en1b = gr[en - 1], en2b = gr[en - 2];

      //if (CanPair(stb, enb)) {  // TODO lonely pairs?
      if (energy::ViableFoldingPair(st, en)) {
        penergy_t p{0};
        const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
        for (int ist = st + 1; ist < st + max_inter + 2; ++ist)  // TODO lyngso's ?
          for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien)
            p += FastTwoLoop(st, en, ist, ien) * gpt[ist][ien][PT_P];
        // Hairpin loops.
        p += FastHairpin(st, en);
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        const penergy_t base_branch_cost = gppc.augubranch[stb][enb] * gppc.em.multiloop_hack_a;

        // (<   ><   >)
        p += base_branch_cost * gpt[st + 1][en - 1][PT_U2];
        // (3<   ><   >) 3'
        p += base_branch_cost * gpt[st + 2][en - 1][PT_U2] * gppc.em.dangle3[stb][st1b][enb];
        // (<   ><   >5) 5'
        p += base_branch_cost * gpt[st + 1][en - 2][PT_U2] * gppc.em.dangle5[stb][en1b][enb];
        // (.<   ><   >.) Terminal mismatch
        p += base_branch_cost * gpt[st + 2][en - 2][PT_U2] *
            gppc.em.terminal[stb][st1b][en1b][enb];

        const auto outer_coax = gppc.em.MismatchCoaxial(stb, st1b, en1b, enb);
        for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
          // Paired coaxial stacking cases:
          base_t pl1b = gr[piv - 1], plb = gr[piv], prb = gr[piv + 1], pr1b = gr[piv + 2];

          // (.(   )   .) Left outer coax - P
          p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[piv + 1][en - 2][PT_U] *
              gppc.augubranch[st2b][plb] * outer_coax;
          // (.   (   ).) Right outer coax
          p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[piv + 1][en - 2][PT_P] *
              gppc.augubranch[prb][en2b] * outer_coax;

          // (.(   ).   ) Left right coax
          p += base_branch_cost * gpt[st + 2][piv - 1][PT_P] * gpt[piv + 1][en - 1][PT_U] *
              gppc.augubranch[st2b][pl1b] * gppc.em.MismatchCoaxial(pl1b, plb, st1b, st2b);
          // (   .(   ).) Right left coax
          p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 2][en - 2][PT_P] *
              gppc.augubranch[pr1b][en2b] * gppc.em.MismatchCoaxial(en2b, en1b, prb, pr1b);

          // ((   )   ) Left flush coax
          p += base_branch_cost * gpt[st + 1][piv][PT_P] * gpt[piv + 1][en - 1][PT_U] *
              gppc.augubranch[st1b][plb] * gppc.em.stack[stb][st1b][plb][enb];
          // (   (   )) Right flush coax
          p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[piv + 1][en - 1][PT_P] *
              gppc.augubranch[prb][en1b] * gppc.em.stack[stb][prb][en1b][enb];
        }

        gpt[st][en][PT_P] = p;
      }
      penergy_t u{0}, u2{0}, rcoax{0}, wc{0}, gu{0};
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
        const penergy_t base00 = gpt[st][piv][PT_P] * gppc.augubranch[stb][pb];
        const penergy_t base01 = gpt[st][piv - 1][PT_P] * gppc.augubranch[stb][pl1b];
        const penergy_t base10 = gpt[st + 1][piv][PT_P] * gppc.augubranch[st1b][pb];
        const penergy_t base11 = gpt[st + 1][piv - 1][PT_P] * gppc.augubranch[st1b][pl1b];

        // (   )<   > - U, U_WC?, U_GU?
        u2 += base00 * gpt[piv + 1][en][PT_U];
        penergy_t val = base00 + base00 * gpt[piv + 1][en][PT_U];
        u += val;
        if (IsGu(stb, pb)) gu += val;
        else wc += val;

        // (   )3<   > 3' - U
        val = base01 * gppc.em.dangle3[pl1b][pb][stb];
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // 5(   )<   > 5' - U
        val = base10 * gppc.em.dangle5[pb][stb][st1b];
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<   > Terminal mismatch - U
        val = base11 * gppc.em.terminal[pl1b][pb][stb][st1b];
        u += val;
        val *= gpt[piv + 1][en][PT_U];
        u += val;
        u2 += val;

        // .(   ).<(   ) > Left coax - U
        val = base11 * gppc.em.MismatchCoaxial(pl1b, pb, stb, st1b);
        val = val * (gpt[piv + 1][en][PT_U_WC] + gpt[piv + 1][en][PT_U_GU]);
        u += val;
        u2 += val;

        // (   )<.(   ). > Right coax forward and backward
        val = base00 * gpt[piv + 1][en][PT_U_RCOAX];
        u += val;
        u2 += val;
        val = base11 * gppc.em.MismatchCoaxial(pl1b, pb, stb, st1b);
        rcoax += val;
        rcoax += val * gpt[piv + 1][en][PT_U];

        // (   )(<   ) > Flush coax - U
        val = base01 * gppc.em.stack[pl1b][pb][pb ^ 3][stb] * gpt[piv][en][PT_U_WC];
        u += val;
        u2 += val;
        if (pb == G || pb == U) {
          val = base01 * gppc.em.stack[pl1b][pb][pb ^ 1][stb] * gpt[piv][en][PT_U_GU];
          u += val;
          u2 += val;
        }
      }
      gpt[st][en][PT_U] = u;
      gpt[st][en][PT_U2] = u2;
      gpt[st][en][PT_U_WC] = wc;
      gpt[st][en][PT_U_GU] = gu;
      gpt[st][en][PT_U_RCOAX] = rcoax;
    }
  }

  // Compute the exterior tables.
  Exterior();

  // Fill the left triangle.
  // The meaning of the tables changes here:
  // U, U2, WC, GU, RCOAX: any table index with en < st must have a loop enclosing (en, st)
  for (int st = N - 1; st >= 0; --st) {
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
        penergy_t p{0};
        const int ost_max = std::min(st + TWOLOOP_MAX_SZ + 2, N);
        for (int ost = st + 1; ost < ost_max; ++ost) {
          const int oen_min = std::max(en - TWOLOOP_MAX_SZ - 1 + (ost - st - 1), 0);
          for (int oen = en - 1; oen >= oen_min; --oen)
            p += FastTwoLoop(oen, ost, en, st) * gpt[ost][oen][PT_P];
        }
        const penergy_t base_branch_cost = gppc.augubranch[stb][enb] * gppc.em.multiloop_hack_a;
        const penergy_t outer_coax = lspace && rspace ?
            gppc.em.MismatchCoaxial(stb, st1b, en1b, enb) : 0.0;
        // Try being an exterior loop - coax cases handled in the loop after this.
        {
          const penergy_t augu = gppc.em.AuGuPenalty(enb, stb);
          const penergy_t rext = gptext[st + 1][PTEXT_R];
          const penergy_t r1ext = lspace > 1 ? gptext[st + 2][PTEXT_R] : penergy_t{1};
          const penergy_t lext = rspace ? gptext[en - 1][PTEXT_L] : penergy_t{1};
          const penergy_t l1ext = rspace > 1 ? gptext[en - 2][PTEXT_L] : penergy_t{1};

          // |<   >)   (<   >| - Exterior loop
          p += augu * lext * rext;
          if (lspace) {
            // |<   >(   )3<   >| 3' - Exterior loop
            // lspace > 0
            p += augu * lext * r1ext * gppc.em.dangle3[stb][st1b][enb];
            // |  >5)   (<   | 5' - Enclosing loop
            if (rspace > 1)
              p += base_branch_cost * gpt[st + 1][en - 2][PT_U2] * gppc.em.dangle5[stb][en1b][enb];
          }
          if (rspace) {
            // |<   >5(   )<   >| 5' - Exterior loop
            // rspace > 0
            p += augu * l1ext * rext * gppc.em.dangle5[stb][en1b][enb];
            // |   >)   (3<  | 3' - Enclosing loop
            if (lspace > 1)
              p += base_branch_cost * gpt[st + 2][en - 1][PT_U2] * gppc.em.dangle3[stb][st1b][enb];
          }
          if (lspace && rspace) {
            // |<   >m(   )m<   >| Terminal mismatch - Exterior loop
            // lspace > 0 && rspace > 0
            p += augu * l1ext * r1ext * gppc.em.terminal[stb][st1b][en1b][enb];
            // |   >)   (<   | - Enclosing loop
            p += base_branch_cost * gpt[st + 1][en - 1][PT_U2];
          }
          // |  >m)   (m<  | Terminal mismatch - Enclosing loop
          if (lspace > 1 && rspace > 1)
            p += base_branch_cost * gpt[st + 2][en - 2][PT_U2] *
                gppc.em.terminal[stb][st1b][en1b][enb];

          const int limit = en + N;
          for (int tpiv = st; tpiv <= limit; ++tpiv) {
            const int pl = FastMod(tpiv - 1, N), piv = FastMod(tpiv, N), pr = FastMod(tpiv + 1, N);
            base_t pl1b = gr[pl], plb = gr[piv], prb = gr[pr];
            const bool left_formable = tpiv < N && tpiv - st - 2 >= HAIRPIN_MIN_SZ;
            const bool right_formable = tpiv >= N && en - piv - 2 >= HAIRPIN_MIN_SZ;

            if (left_formable) {
              const penergy_t rp1ext = gptext[piv + 1][PTEXT_R];
              // |<   >)   (.(   ).<   >| Exterior loop - Left right coax
              // lspace > 1 && not enclosed
              p += gpt[st + 2][pl][PT_P] * lext * rp1ext * gppc.em.AuGuPenalty(stb, enb) *
                  gppc.em.AuGuPenalty(st2b, pl1b) * gppc.em.MismatchCoaxial(pl1b, plb, st1b, st2b);
              // |<   >)   ((   )<   >| Exterior loop - Left flush coax
              // lspace > 0 && not enclosed
              p += gpt[st + 1][piv][PT_P] * lext * rp1ext * gppc.em.AuGuPenalty(stb, enb) *
                  gppc.em.AuGuPenalty(st1b, plb) * gppc.em.stack[stb][st1b][plb][enb];

              if (rspace) {
                // |<   >.)   (.(   )<   >| Exterior loop - Left outer coax
                // lspace > 0 && rspace > 0 && not enclosed
                p += gpt[st + 2][piv][PT_P] * l1ext * rp1ext * gppc.em.AuGuPenalty(stb, enb) *
                    gppc.em.AuGuPenalty(st2b, plb) * outer_coax;
              }
            }

            if (right_formable) {
              const penergy_t lpext = piv > 0 ? gptext[piv - 1][PTEXT_L] : penergy_t{1};
              // |<   >.(   ).)   (<   >| Exterior loop - Right left coax
              // rspace > 1 && not enclosed
              p += gpt[pr][en - 2][PT_P] * lpext * rext * gppc.em.AuGuPenalty(stb, enb) *
                  gppc.em.AuGuPenalty(prb, en2b) * gppc.em.MismatchCoaxial(en2b, en1b, plb, prb);
              // |<   >(   ))   (<   >| Exterior loop - Right flush coax
              // rspace > 0 && not enclosed
              p += gpt[piv][en - 1][PT_P] * lpext * rext * gppc.em.AuGuPenalty(stb, enb) *
                  gppc.em.AuGuPenalty(plb, en1b) * gppc.em.stack[en1b][enb][stb][plb];

              if (lspace) {
                // |<   >(   ).)   (.<   >| Exterior loop - Right outer coax
                // lspace > 0 && rspace > 1 && not enclosed
                p += gpt[piv][en - 2][PT_P] * lpext * r1ext * gppc.em.AuGuPenalty(stb, enb) *
                    gppc.em.AuGuPenalty(plb, en2b) * outer_coax;
              }
            }
          }
        }

        // Enclosing loop cases.
        // Can start at st + 2 because we need to forman enclosing loop.
        // At worst the enclosing loop has to start at st + 1.
        const int limit = en + N;
        for (int tpiv = st + 2; tpiv < limit; ++tpiv) {
          const int pl = FastMod(tpiv - 1, N), piv = FastMod(tpiv, N),
              pr = FastMod(tpiv + 1, N), pr1 = FastMod(tpiv + 2, N);
          // Left block is: [st, piv], Right block is: [piv + 1, en].
          base_t pl1b = gr[pl], plb = gr[piv], prb = gr[pr], pr1b = gr[pr1];
          // When neither the left block nor the right block straddles the border
          // we don't get enclosing loops sometimes, so check against N - 1.
          const bool straddling = tpiv != (N - 1);
          // Left loop formable if not straddling, and is big enough or crosses over.
          const bool left_formable = straddling && (piv - st - 2 >= HAIRPIN_MIN_SZ || tpiv >= N);
          const bool right_formable = straddling && (en - piv - 3 >= HAIRPIN_MIN_SZ || tpiv < N - 1);
          const bool left_dot_formable = left_formable && tpiv != N;  // Can't split a dot.
          const bool right_dot_formable = right_formable && tpiv != N - 2;  // Can't split a dot.

          if (lspace > 1 && rspace > 1) {
            if (left_formable) {
              // |  >.)   (.(   )<  | Enclosing loop - Left outer coax
              // lspace > 1 && rspace > 1 && enclosed
              p += base_branch_cost * gpt[st + 2][piv][PT_P] * gpt[pr][en - 2][PT_U] *
                  gppc.augubranch[st2b][plb] * outer_coax;
            }

            if (right_formable) {
              // |  >(   ).)   (.<  | Enclosing loop - Right outer coax
              // lspace > 1 && rspace > 1 && enclosed
              p += base_branch_cost * gpt[st + 2][piv][PT_U] * gpt[pr][en - 2][PT_P] *
                  gppc.augubranch[prb][en2b] * outer_coax;
            }
          }

          if (lspace > 1 && rspace && left_dot_formable) {
            // |  >)   (.(   ).<  | Enclosing loop - Left right coax
            // lspace > 1 && rspace > 0 && enclosed && no dot split
            p += base_branch_cost * gpt[st + 2][pl][PT_P] * gpt[pr][en - 1][PT_U] *
                gppc.augubranch[st2b][pl1b] * gppc.em.MismatchCoaxial(pl1b, plb, st1b, st2b);
          }

          if (lspace && rspace > 1 && right_dot_formable) {
            // |  >.(   ).)   (<  | Enclosing loop - Right left coax
            // lspace > 0 && rspace > 1 && enclosed && no dot split
            p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[pr1][en - 2][PT_P] *
                gppc.augubranch[pr1b][en2b] * gppc.em.MismatchCoaxial(en2b, en1b, prb, pr1b);
          }

          if (lspace && rspace) {
            if (left_formable) {
              // |  >)   ((   )<  | Enclosing loop - Left flush coax
              // lspace > 0 && rspace > 0 && enclosed
              p += base_branch_cost * gpt[st + 1][piv][PT_P] * gpt[pr][en - 1][PT_U] *
                  gppc.augubranch[st1b][plb] * gppc.em.stack[stb][st1b][plb][enb];
            }

            if (right_formable) {
              // |  >(   ))   (<  | Enclosing loop - Right flush coax
              // lspace > 0 && rspace > 0 && enclosed
              p += base_branch_cost * gpt[st + 1][piv][PT_U] * gpt[pr][en - 1][PT_P] *
                  gppc.augubranch[prb][en1b] * gppc.em.stack[en1b][enb][stb][prb];
            }
          }
        }
        gpt[st][en][PT_P] = p;
      }
      penergy_t u{0}, u2{0}, rcoax{0}, wc{0}, gu{0};
      // Update unpaired.
      // Choose |st| to be unpaired, but only if we can maintain the constraint that we have
      // an enclosing loop formed.
      if (st + 1 < N) {
        u += gpt[st + 1][en][PT_U];
        u2 += gpt[st + 1][en][PT_U2];
      }

      for (int tpiv = st + 1; tpiv <= en + N; ++tpiv) {
        const int pl = FastMod(tpiv - 1, N), piv = FastMod(tpiv, N), pr = FastMod(tpiv + 1, N);
        const auto pb = gr[piv], pl1b = gr[pl], prb = gr[pr];
        const penergy_t base00 = gpt[st][piv][PT_P] * gppc.augubranch[stb][pb];
        const penergy_t base01 = gpt[st][pl][PT_P] * gppc.augubranch[stb][pl1b];

        // Must have an enclosing loop.
        const bool straddling = tpiv != N - 1;
        const bool dot_straddling = straddling && tpiv != N;
        penergy_t val{0};

        if (straddling) {
          // |  >>   <(   )<  |
          val = base00 * gpt[pr][en][PT_U];
          u2 += val;
          u += val;
          if (IsGu(stb, pb)) gu += val;
          else wc += val;
          // U must cross the boundary to have the rest of it be nothing.
          if (tpiv >= N) {
            u += base00;
            if (IsGu(stb, pb)) gu += base00;
            else wc += base00;
          }

          // |  )  >>   <(   )<(  | Flush coax
          // straddling
          val = base00 * gppc.em.stack[pb][prb][prb ^ 3][stb] * gpt[pr][en][PT_U_WC];
          u += val;
          u2 += val;
          if (prb == G || prb == U) {
            val = base00 * gppc.em.stack[pb][prb][prb ^ 1][stb] * gpt[pr][en][PT_U_GU];
            u += val;
            u2 += val;
          }
          // |  ).  >>   <(   )<.(   | Right coax forward
          // straddling
          val = base00 * gpt[pr][en][PT_U_RCOAX];
          u += val;
          u2 += val;
        }

        if (dot_straddling) {
          // |  >>   <(   )3<  | 3'
          val = base01 * gppc.em.dangle3[pl1b][pb][stb];
          if (tpiv >= N) u += val;
          val *= gpt[pr][en][PT_U];
          u += val;
          u2 += val;
        }

        if (lspace) {
          const penergy_t base10 = gpt[st + 1][piv][PT_P] * gppc.augubranch[st1b][pb];
          const penergy_t base11 = gpt[st + 1][pl][PT_P] * gppc.augubranch[st1b][pl1b];

          if (straddling) {
            // |  >>   <5(   )<  | 5'
            val = base10 * gppc.em.dangle5[pb][stb][st1b];
            if (tpiv >= N) u += val;
            val *= gpt[pr][en][PT_U];
            u += val;
            u2 += val;
          }

          if (dot_straddling) {
            // |  >>   <.(   ).<  | Terminal mismatch
            // lspace > 0 && dot_straddling
            val = base11 * gppc.em.terminal[pl1b][pb][stb][st1b];
            if (tpiv >= N) u += val;
            val *= gpt[pr][en][PT_U];
            u += val;
            u2 += val;
            // |  )>>   <.(   ).<(  | Left coax
            // lspace > 0 && dot_straddling
            val = base11 * gppc.em.MismatchCoaxial(pl1b, pb, stb, st1b);
            val = val * (gpt[pr][en][PT_U_WC] + gpt[pr][en][PT_U_GU]);
            u += val;
            u2 += val;
            // |  ).  >>   <(   )<.(   | Right coax backward
            val = base11 * gppc.em.MismatchCoaxial(pl1b, pb, stb, st1b);
            if (tpiv >= N) rcoax += val;
            rcoax += val * gpt[pr][en][PT_U];
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