#!/usr/bin/env python3
# Copyright 2016, Eliot Courtney.
#
# This file is part of kekrna.
#
# kekrna is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# kekrna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with kekrna.
# If not, see <http://www.gnu.org/licenses/>.
from common import *

icam1 = 'GCGCCCCAGUCGACGCUGAGCUCCUCUGCUACUCAGAGUUGCAACCUCAGCCUCGCUAUGGCUCCCAGCAGCCCCCGGCCCGCGCUGCCCGCACUCCUGGUCCUGCUCGGGGCUCUGUUCCCAGGACCUGGCAAUGCCCAGACAUCUGUGUCCCCCUCAAAAGUCAUCCUGCCCCGGGGAGGCUCCGUGCUGGUGACAUGCAGCACCUCCUGUGACCAGCCCAAGUUGUUGGGCAUAGAGACCCCGUUGCCUAAAAAGGAGUUGCUCCUGCCUGGGAACAACCGGAAGGUGUAUGAACUGAGCAAUGUGCAAGAAGAUAGCCAACCAAUGUGCUAUUCAAACUGCCCUGAUGGGCAGUCAACAGCUAAAACCUUCCUCACCGUGUACUGGACUCCAGAACGGGUGGAACUGGCACCCCUCCCCUCUUGGCAGCCAGUGGGCAAGAACCUUACCCUACGCUGCCAGGUGGAGGGUGGGGCACCCCGGGCCAACCUCACCGUGGUGCUGCUCCGUGGGGAGAAGGAGCUGAAACGGGAGCCAGCUGUGGGGGAGCCCGCUGAGGUCACGACCACGGUGCUGGUGAGGAGAGAUCACCAUGGAGCCAAUUUCUCGUGCCGCACUGAACUGGACCUGCGGCCCCAAGGGCUGGAGCUGUUUGAGAACACCUCGGCCCCCUACCAGCUCCAGACCUUUGUCCUGCCAGCGACUCCCCCACAACUUGUCAGCCCCCGGGUCCUAGAGGUGGACACGCAGGGGACCGUGGUCUGUUCCCUGGACGGGCUGUUCCCAGUCUCGGAGGCCCAGGUCCACCUGGCACUGGGGGACCAGAGGUUGAACCCCACAGUCACCUAUGGCAACGACUCCUUCUCGGCCAAGGCCUCAGUCAGUGUGACCGCAGAGGACGAGGGCACCCAGCGGCUGACGUGUGCAGUAAUACUGGGGAACCAGAGCCAGGAGACACUGCAGACAGUGACCAUCUACAGCUUUCCGGCGCCCAACGUGAUUCUGACGAAGCCAGAGGUCUCAGAAGGGACCGAGGUGACAGUGAAGUGUGAGGCCCACCCUAGAGCCAAGGUGACGCUGAAUGGGGUUCCAGCCCAGCCACUGGGCCCGAGGGCCCAGCUCCUGCUGAAGGCCACCCCAGAGGACAACGGGCGCAGCUUCUCCUGCUCUGCAACCCUGGAGGUGGCCGGCCAGCUUAUACACAAGAACCAGACCCGGGAGCUUCGUGUCCUGUAUGGCCCCCGACUGGACGAGAGGGAUUGUCCGGGAAACUGGACGUGGCCAGAAAAUUCCCAGCAGACUCCAAUGUGCCAGGCUUGGGGGAACCCAUUGCCCGAGCUCAAGUGUCUAAAGGAUGGCACUUUCCCACUGCCCAUCGGGGAAUCAGUGACUGUCACUCGAGAUCUUGAGGGCACCUACCUCUGUCGGGCCAGGAGCACUCAAGGGGAGGUCACCCGCGAGGUGACCGUGAAUGUGCUCUCCCCCCGGUAUGAGAUUGUCAUCAUCACUGUGGUAGCAGCCGCAGUCAUAAUGGGCACUGCAGGCCUCAGCACGUACCUCUAUAACCGCCAGCGGAAGAUCAAGAAAUACAGACUACAACAGGCCCAAAAAGGGACCCCCAUGAAACCGAACACACAAGCCACGCCUCCCUGAACCUAUCCCGGGACAGGGCCUCUUCCUCGGCCUUCCCAUAUUGGUGGCAGUGGUGCCACACUGAACAGAGUGGAAGACAUAUGCCAUGCAGCUACACCUACCGGCCCUGGGACGCCGGAGGACAGGGCAUUGUCCUCAGUCAGAUACAACAGCAUUUGGGGCCAUGGUACCUGCACACCUAAAACACUAGGCCACGCAUCUGAUCUGUAGUCACAUGACUAAGCCAAGAGGAAGGAGCAAGACUCAAGACAUGAUUGAUGGAUGUUAAAGUCUAGCCUGAUGAGAGGGGAAGUGGUGGGGGAGACAUAGCCCCACCAUGAGGACAUACAACUGGGAAAUACUGAAACUUGCUGCCUAUUGGGUAUGCUGAGGCCCACAGACUUACAGAAGAAGUGGCCCUCCAUAGACAUGUGUAGCAUCAAAACACAAAGGCCCACACUUCCUGACGGAUGCCAGCUUGGGCACUGCUGUCUACUGACCCCAACCCUUGAUGAUAUGUAUUUAUUCAUUUGUUAUUUUACCAGCUAUUUAUUGAGUGUCUUUUAUGUAGGCUAAAUGAACAUAGGUCUCUGGCCUCACGGAGCUCCCAGUCCAUGUCACAUUCAAGGUCACCAGGUACAGUUGUACAGGUUGUACACUGCAGGAGAGUGCCUGGCAAAAAGAUCAAAUGGGGCUGGGACUUCUCAUUGGCCAACCUGCCUUUCCCCAGAAGGAGUGAUUUUUCUAUCGGCACAAAAGCACUAUAUGGACUGGUAAUGGUUCACAGGUUCAGAGAUUACCCAGUGAGGCCUUAUUCCUCCCUUCCCCCCAAAACUGACACCUUUGUUAGCCACCUCCCCACCCACAUACAUUUCUGCCAGUGUUCACAAUGACACUCAGCGGUCAUGUCUGGACAUGAGUGCCCAGGGAAUAUGCCCAAGCUAUGCCUUGUCCUCUUGUCCUGUUUGCAUUUCACUGGGAGCUUGCACUAUUGCAGCUCCAGUUUCCUGCAGUGAUCAGGGUCCUGCAAGCAGUGGGGAAGGGGGCCAAGGUAUUGGAGGACUCCCUCCCAGCUUUGGAAGGGUCAUCCGCGUGUGUGUGUGUGUGUAUGUGUAGACAAGCUCUCGCUCUGUCACCCAGGCUGGAGUGCAGUGGUGCAAUCAUGGUUCACUGCAGUCUUGACCUUUUGGGCUCAAGUGAUCCUCCCACCUCAGCCUCCUGAGUAGCUGGGACCAUAGGCUCACAACACCACACCUGGCAAAUUUGAUUUUUUUUUUUUUUUUCAGAGACGGGGUCUCGCAACAUUGCCCAGACUUCCUUUGUGUUAGUUAAUAAAGCUUUCUCAACUGCC'
lendeltanum = [
  (1000, 6, 500000), (1200, 6, 500000), (1600, 6, 500000), (1700, 6, 4000000)]


def run_num(alg, ldn):
  print('alg%s - nums, quiet' % alg)
  for l, _, num in ldn:
    res = run_command(
      os.path.join('build', 'c++-release', 'subopt'),
      '-q', '-num', str(num), '-subopt-alg', alg, icam1[:l])
    print('  len %d, num %d: %s' % (l, num, res))


def run_delta(alg, ldn):
  print('alg%s - deltas, quiet' % alg)
  for l, d, _ in ldn:
    res = run_command(
      os.path.join('build', 'c++-release', 'subopt'),
      '-q', '-delta', str(d), '-subopt-alg', alg, icam1[:l])
    print('  len %d, delta %d: %s' % (l, d, res))


def run_fold(alg, lens):
  print('alg%s -folding' % alg)
  for l in lens:
    res = run_command(
      os.path.join('build', 'c++-release', 'fold'),
      '-dp-alg', alg, icam1[:l])
    print('  len %d: %s' % (l, res))

def main():
  run_fold('2', [1800, 2500, 2900])
  # run_delta('1', lendeltanum)
  # run_num('1', lendeltanum)


if __name__ == '__main__':
  main()
