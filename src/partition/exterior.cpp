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
#include "partition/partition_globals.h"
#include "energy/energy_globals.h"
#include "energy/fast_energy.h"

namespace kekrna {
namespace partition {
namespace internal {

using energy::gem;
using energy::Boltzmann;

void Exterior() {
  const int N = int(gr.size());
  gptext[N][PTEXT_R] = 1.0;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    gptext[st][PTEXT_R] += gptext[st + 1][PTEXT_R];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const auto base00 = gpt[st][en][PT_P] * Boltzmann(gem.AuGuPenalty(stb, enb));
      const auto base01 = gpt[st][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(stb, en1b));
      const auto base10 = gpt[st + 1][en][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, enb));
      const auto base11 = gpt[st + 1][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, en1b));

      // (   )<   >
      auto val = base00 * gptext[en + 1][PTEXT_R];
      gptext[st][PTEXT_R] += val;
      if (IsGu(stb, enb)) gptext[st][PTEXT_R_GU] += val;
      else gptext[st][PTEXT_R_WC] += val;

      // (   )3<   > 3'
      gptext[st][PTEXT_R] += base01 * Boltzmann(gem.dangle3[en1b][enb][stb]) * gptext[en + 1][PTEXT_R];
      // 5(   )<   > 5'
      gptext[st][PTEXT_R] += base10 * Boltzmann(gem.dangle5[enb][stb][st1b]) * gptext[en + 1][PTEXT_R];
      // .(   ).<   > Terminal mismatch
      gptext[st][PTEXT_R] += base11 *
          Boltzmann(gem.terminal[en1b][enb][stb][st1b]) * gptext[en + 1][PTEXT_R];
      // .(   ).<(   ) > Left coax
      val = base11 * Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b));
      gptext[st][PTEXT_R] += val * gptext[en + 1][PTEXT_R_GU];
      gptext[st][PTEXT_R] += val * gptext[en + 1][PTEXT_R_WC];

      // (   ).<(   ). > Right coax forward
      gptext[st][PTEXT_R] += base01 * gptext[en + 1][PTEXT_R_RCOAX];
      // (   ).<( * ). > Right coax backward
      if (st > 0)
        gptext[st][PTEXT_R_RCOAX] += base01 *
            Boltzmann(gem.MismatchCoaxial(en1b, enb, gr[st - 1], stb)) * gptext[en + 1][PTEXT_R];

      if (en < N - 1) {
        // (   )<(   ) > Flush coax
        const auto enrb = gr[en + 1];
        gptext[st][PTEXT_R] += base00 *
            Boltzmann(gem.stack[enb][enrb][enrb ^ 3][stb]) * gptext[en + 1][PTEXT_R_WC];
        if (enrb == G || enrb == U)
          gptext[st][PTEXT_R] += base00 *
              Boltzmann(gem.stack[enb][enrb][enrb ^ 1][stb]) * gptext[en + 1][PTEXT_R_GU];
      }
    }
  }

  gptext[0][PTEXT_L] = 1.0;
  for (int en = 1; en < N; ++en) {
    // Case: No pair ending here
    gptext[en][PTEXT_L] += gptext[en - 1][PTEXT_L];
    for (int st = 0; st < en - HAIRPIN_MIN_SZ; ++st) {
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const auto base00 = gpt[st][en][PT_P] * Boltzmann(gem.AuGuPenalty(stb, enb));
      const auto base01 = gpt[st][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(stb, en1b));
      const auto base10 = gpt[st + 1][en][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, enb));
      const auto base11 = gpt[st + 1][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, en1b));
      const auto ptextl = st ? gptext[st - 1][PTEXT_L] : 1.0;
      const auto ptextlgu = st ? gptext[st - 1][PTEXT_L_GU] : 0.0;
      const auto ptextlwc = st ? gptext[st - 1][PTEXT_L_WC] : 0.0;
      const auto ptextllcoaxx = st ? gptext[st - 1][PTEXT_L_LCOAX] : 0.0;

      // <   >(   )
      auto val = base00 * ptextl;
      gptext[en][PTEXT_L] += val;
      if (IsGu(stb, enb)) gptext[en][PTEXT_L_GU] += val;
      else gptext[en][PTEXT_L_WC] += val;

      // <   >(   )3 3'
      gptext[en][PTEXT_L] += base01 * Boltzmann(gem.dangle3[en1b][enb][stb]) * ptextl;
      // 5(   )<   > 5'
      gptext[en][PTEXT_L] += base10 * Boltzmann(gem.dangle5[enb][stb][st1b]) * ptextl;
      // .(   ).<   > Terminal mismatch
      gptext[en][PTEXT_L] += base11 * Boltzmann(gem.terminal[en1b][enb][stb][st1b]) * ptextl;
      // <  (   )>.(   ). Right coax
      val = base11 * Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b));
      gptext[en][PTEXT_L] += val * ptextlgu;
      gptext[en][PTEXT_L] += val * ptextlwc;

      // <  .(   ).>(   ) Left coax forward
      gptext[en][PTEXT_L] += base00 * ptextllcoaxx;
      // <  .( * ).>(   ) Left coax backward
      gptext[en][PTEXT_L_LCOAX] += base11 *
          Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b)) * ptextl;

      if (st) {
        // < (   )>(   ) Flush coax
        const auto stl1b = gr[st - 1];
        gptext[en][PTEXT_L] += base00 *
            Boltzmann(gem.stack[stl1b][stb][stl1b ^ 3][enb]) * ptextlwc;
        if (stl1b == G || stl1b == U)
          gptext[en][PTEXT_L] += base00 *
              Boltzmann(gem.stack[stl1b][stb][stl1b ^ 1][enb]) * ptextlgu;
      }
    }
  }

  assert(std::abs(gptext[N - 1][PTEXT_L] - gptext[0][PTEXT_R]) < 1e-6);
}

}
}
}
