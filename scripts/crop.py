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


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('files', nargs='+')
  args = parser.parse_args()
  for i in args.files:
    os.system('convert %s -trim %s' % (i, i))


if __name__ == '__main__':
  main()
