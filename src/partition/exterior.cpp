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
  gptext[N][PTEXT] = 1.0;
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    gptext[st][PTEXT] += gptext[st + 1][PTEXT];
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      const auto stb = gr[st], st1b = gr[st + 1], enb = gr[en], en1b = gr[en - 1];
      const auto base00 = gpt[st][en][PT_P] * Boltzmann(gem.AuGuPenalty(stb, enb));
      const auto base01 = gpt[st][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(stb, en1b));
      const auto base10 = gpt[st + 1][en][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, enb));
      const auto base11 = gpt[st + 1][en - 1][PT_P] * Boltzmann(gem.AuGuPenalty(st1b, en1b));

      // (   )<   >
      auto val = base00 * gptext[en + 1][PTEXT];
      gptext[st][PTEXT] += val;
      if (IsGu(stb, enb)) gptext[st][PTEXT_GU] += val;
      else gptext[st][PTEXT_WC] += val;

      // (   )3<   > 3'
      gptext[st][PTEXT] += base01 * Boltzmann(gem.dangle3[en1b][enb][stb]) * gptext[en + 1][PTEXT];
      // 5(   )<   > 5'
      gptext[st][PTEXT] += base10 * Boltzmann(gem.dangle5[enb][stb][st1b]) * gptext[en + 1][PTEXT];
      // .(   ).<   > Terminal mismatch
      gptext[st][PTEXT] += base11 *
          Boltzmann(gem.terminal[en1b][enb][stb][st1b]) * gptext[en + 1][PTEXT];
      // .(   ).<(   ) > Left coax  x
      val = base11 * Boltzmann(gem.MismatchCoaxial(en1b, enb, stb, st1b));
      gptext[st][PTEXT] += val * gptext[en + 1][PTEXT_GU];
      gptext[st][PTEXT] += val * gptext[en + 1][PTEXT_WC];

      // (   ).<(   ). > Right coax forward
      gptext[st][PTEXT] += base01 * gptext[en + 1][PTEXT_RCOAX];
      // (   ).<( * ). > Right coax backward
      if (st > 0)
        gptext[st][PTEXT_RCOAX] += base01 *
            Boltzmann(gem.MismatchCoaxial(en1b, enb, gr[st - 1], stb)) * gptext[en + 1][PTEXT];

      if (en < N - 1) {
        // (   )<(   ) > Flush coax
        const auto enrb = gr[en + 1];
        gptext[st][PTEXT] += base00 *
            Boltzmann(gem.stack[enb][enrb][enrb ^ 3][stb]) * gptext[en + 1][PTEXT_WC];
        if (enrb == G || enrb == U)
          gptext[st][PTEXT] += base00 *
              Boltzmann(gem.stack[enb][enrb][enrb ^ 1][stb]) * gptext[en + 1][PTEXT_GU];
      }
    }
  }
}

}
}
}
