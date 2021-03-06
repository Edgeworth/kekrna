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
import argparse
import os
import sys
import shutil

parser = argparse.ArgumentParser()
parser.add_argument(
  '-t', '--type', choices=['debug', 'asan', 'ubsan', 'release', 'relwithdebinfo'],
  default='debug', required=False)
parser.add_argument('-c', '--use_clang', action='store_true', default=False, required=False)
parser.add_argument('-a', '--use_afl', action='store_true', default=False, required=False)
parser.add_argument('-d', '--dry', action='store_true', default=False, required=False)
parser.add_argument('--compilers', type=str, nargs=2, required=False)
parser.add_argument('-r', '--regenerate', action='store_true', default=False, required=False)
parser.add_argument('targets', nargs='*', type=str)
args = parser.parse_args()

def run_command(cmd):
  if args.dry:
    print(cmd)
  else:
    res = os.system(cmd)
    if res != 0:
      sys.exit(1)

if bool(args.use_clang) + bool(args.use_afl) + bool(args.compilers) > 1:
  parser.error('At most one compiler related flag allowed simultaneously')

compilers = ('cc', 'c++')

if args.use_clang:
  compilers = ('clang', 'clang++')
elif args.use_afl:
  compilers = ('afl-clang-fast', 'afl-clang-fast++')
elif args.compilers:
  compilers = args.compilers

defs = {
  'CMAKE_C_COMPILER': compilers[0],
  'CMAKE_CXX_COMPILER': compilers[1],
  'CMAKE_BUILD_TYPE': args.type
}

build_dir = os.path.join('build', defs['CMAKE_CXX_COMPILER'] + '-' + defs['CMAKE_BUILD_TYPE'])
regenerate = args.regenerate

if regenerate and os.path.exists(build_dir):
  shutil.rmtree(build_dir)

if not os.path.exists(build_dir):
  os.makedirs(build_dir)
  regenerate = True

os.chdir(build_dir)
if regenerate:
  print('Regenerating cmake files.')
  def_str = ' '.join('-D %s=\'%s\'' % (i, k) for i, k in defs.items())
  run_command('cmake %s ../../' % def_str)
run_command('AFL_HARDEN=1 make -j16 %s' % ' '.join(args.targets))
os.chdir('../../')
