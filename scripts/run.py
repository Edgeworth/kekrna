#!/usr/bin/env python3
import argparse
import shutil
import tempfile

from common import *
from kekvault import KekVault
from rna import *


class RNAstructureDistribution:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.RNASTRUCTURE_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)
    os.putenv('DATAPATH', os.path.join(self.loc, 'data_tables'))

  def fold(self, rna):
    with tempfile.NamedTemporaryFile('w') as f, tempfile.NamedTemporaryFile('r') as out:
      prev_dir = os.getcwd()
      os.chdir(self.loc)
      f.write(rna.to_seq_file())
      f.flush()
      benchmark_results, _ = benchmark_command(
        os.path.join(self.loc, 'exe', 'Fold'), '-mfe', f.name, out.name)
      output = out.read()
      predicted = RNA.from_any_file(output)
      os.chdir(prev_dir)
    return predicted, benchmark_results

  def efn(self, rna):
    with tempfile.NamedTemporaryFile('w') as f, tempfile.NamedTemporaryFile('r') as out:
      f.write(rna.to_ct_file())
      f.flush()
      # Note that not giving the -s flag doesn't make it logarithmic.
      # RNAstructure 5.8 adds the logarithmic and asymmetry models together in this case.
      # RNAstructure also uses a coefficient of -6 for the number of branches, rather than
      # the fitted -9.
      benchmark_results, _ = benchmark_command(
        os.path.join(self.loc, 'exe', 'efn2'), '-s', f.name, out.name)
      output = out.read()
      match = re.search(r'[eE]nergy = (.+)', output.strip())
      energy = float(match.group(1))
    return energy, benchmark_results

  def suboptimal(self, rna, delta):
    # TODO implement
    raise NotImplementedError

  def close(self):
    pass

  def __str__(self):
    return 'RNAstructureDistribution'


class HarnessFolder:
  def __init__(self, loc, name, flag):
    try:
      import default_paths
      loc = loc or default_paths.KEKRNA_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)
    self.name = name
    self.flag = flag

  def fold(self, rna):
    prev_dir = os.getcwd()
    os.chdir(self.loc)
    benchmark_results, stdout = benchmark_command(
      os.path.join('build', 'c++-release', 'harness'), '-f', self.flag, input=rna.seq)
    os.chdir(prev_dir)
    lines = stdout.strip().split('\n')
    assert len(lines) == 2
    return RNA.from_name_seq_db(rna.name, rna.seq, lines[1]), benchmark_results

  def batch_efn(self, rnas):
    prev_dir = os.getcwd()
    os.chdir(self.loc)
    input = '\n'.join('%s\n%s' % (rna.seq, rna.db()) for rna in rnas)
    benchmark_results, stdout = benchmark_command(
      os.path.join('build', 'c++-release', 'harness'), '-e', self.flag, input=input)
    os.chdir(prev_dir)
    energies = [float(i) / 10.0 for i in stdout.strip().split('\n')]
    return energies, benchmark_results

  def efn(self, rna):
    energies, benchmark_results = self.batch_efn([rna])
    return energies[0], benchmark_results

  def suboptimal(self, rna, delta):
    # TODO implement this
    pass

  def close(self):
    pass

  def __str__(self):
    return self.name


class RNAstructureHarness(HarnessFolder):
  def __init__(self, loc=None):
    super().__init__(loc, 'RNAstructureHarness', '-r')


class Rnark(HarnessFolder):
  def __init__(self, loc=None):
    super().__init__(loc, 'Rnark', '-m')


class KekRNA(HarnessFolder):
  def __init__(self, loc=None):
    super().__init__(loc, 'KekRNA', '-k')


class ViennaRNA:
  def __init__(self, loc=None, d3=False):
    try:
      import default_paths
      loc = loc or default_paths.VIENNARNA_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)
    self.extra_args = ['-d2']
    self.name = 'ViennaRNA-d2'
    if d3:
      self.extra_args = ['-d3']
      self.name = 'ViennaRNA-d3'

  def fold(self, rna):
    with tempfile.NamedTemporaryFile('w') as f:
      f.write(rna.seq)
      f.flush()
      benchmark_results, stdout = benchmark_command(
        os.path.join(self.loc, 'src', 'bin', 'RNAfold'),
        *self.extra_args, '--noPS', '-i', f.name)
      seq, db = stdout.strip().split('\n')
      db = db.split(' ')[0]
      predicted = RNA.from_name_seq_db(rna.name, seq.strip(), db.strip())
    return predicted, benchmark_results

  def efn(self, rna):
    raise NotImplementedError

  def suboptimal(self, rna, delta):
    benchmark_results, stdout = benchmark_command(
      os.path.join(self.loc, 'src', 'bin', 'RNAsubopt'),
      *self.extra_args, '-e', '%.1f' % (delta / 10.0), input=rna.seq)
    print(len(stdout.splitlines()))
    return benchmark_results

  def close(self):
    pass

  def __str__(self):
    return self.name


class UNAFold:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.UNAFOLD_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)
    self.tempdir = tempfile.mkdtemp()
    os.putenv('UNAFOLDDAT', os.path.join(self.loc, 'data'))

  def fold(self, rna):
    prev_dir = os.getcwd()
    os.chdir(self.tempdir)
    with open(os.path.join(self.tempdir, 'rna.seq'), 'w') as f:
      f.write(rna.to_db_file())
      f.flush()
      benchmark_results, _ = benchmark_command(
        os.path.join(self.loc, 'src', 'hybrid-ss-min'), f.name)
      predicted = RNA.from_any_file(read_file(os.path.splitext(f.name)[0] + '.ct'))
    os.chdir(prev_dir)
    return predicted, benchmark_results

  def efn(self, rna):
    prev_dir = os.getcwd()
    os.chdir(self.tempdir)
    with open(os.path.join(self.tempdir, 'rna.seq'), 'w') as f:
      f.write(rna.to_ct_file())
      f.flush()
      benchmark_results, stdout = benchmark_command(
        os.path.join(self.loc, 'src', 'ct-energy'), f.name)
      energy = float(stdout.strip())
    os.chdir(prev_dir)
    return energy, benchmark_results

  def suboptimal(self, rna, delta):
    # TODO implement
    raise NotImplementedError

  def close(self):
    shutil.rmtree(self.tempdir)

  def __str__(self):
    return 'UNAFold'


# TODO implement this
class SparseMFEFold:
  def __init__(self, loc=None):
    try:
      import default_paths
      loc = loc or default_paths.SPARSEMFEFOLD_PATH
    except ImportError:
      pass
    assert loc
    self.loc = fix_path(loc)

  def fold(self, rna):
    benchmark_results, stdout = benchmark_command(
      os.path.join(self.loc, 'src', 'SparseMFEFold'), input=rna.seq)
    seq, db = stdout.strip().split('\n')
    db = db.split(' ')[0]
    predicted = RNA.from_name_seq_db(rna.name, seq.strip(), db.strip())
    return predicted, benchmark_results

  def efn(self, rna):
    raise NotImplementedError

  def suboptimal(self, rna, delta):
    raise NotImplementedError

  def close(self):
    pass

  def __str__(self):
    return 'SparseMFEFold'


def run_fold(program, rna):
  frna, benchmark_results = program.fold(rna)
  print('Folding %s with %s: %s\n  %s' % (
    rna.name, program, frna.db(), benchmark_results))


def run_efn(program, rna):
  energy, benchmark_results = program.efn(rna)
  print('Energy of %s with %s: %f\n  %s' % (
    rna.name, program, energy, benchmark_results))


def run_suboptimal(program, rna, delta):
  pass  # TODO - also suboptimals should return # of thingies


BENCHMARK_NUM_TRIES = 5

def run_benchmark(program, dataset, rnastructure_harness):
  if rnastructure_harness is None:
    print('Need RNAstructure for benchmark efn checking')
    sys.exit(1)
  kekvault = KekVault(dataset)
  print('Benchmarking %s on %s' % (program, dataset))
  with open('%s_%s.results' % (program, dataset), 'w') as f:
    idx = 1
    for rna in kekvault:
      print('Running %s on #%d %s' % (program, idx, rna.name))
      prs = []
      for i in range(BENCHMARK_NUM_TRIES):
        prs.append(program.fold(rna))

      accuracy = RNAAccuracy.from_rna(rna, prs[0][0])
      energy, _ = rnastructure_harness.efn(prs[0][0])
      for i, pr in enumerate(prs):
        predicted, result = pr
        f.write('%s %d %d %.5f %.5f %.5f %.5f %.5f %.5f %.2f\n' % (
          rna.name, i, len(rna.seq), result.real, result.usersys, result.maxrss,
          accuracy.fscore, accuracy.ppv, accuracy.sensitivity, energy
        ))
      idx += 1


def parse_rna_from_args(parser, args):
  if bool(args.path) + bool(args.kekvault) + bool(args.cmd) != 1 and not args.benchmark:
    parser.error('Exactly one primary/secondary sequence required')

  if args.path:
    return RNA.from_any_file(read_file(args.path))
  elif args.kekvault:
    return KekVault('archiveii')[args.kekvault]
  elif args.cmd:
    if args.fold:
      if len(args.cmd) != 1:
        parser.error('Direct specification requires one argument for prediction.')
      seq = args.cmd[0]
      return RNA('cmdline', seq, [-1] * len(seq))
    elif args.energy:
      if len(args.cmd) != 2:
        parser.error('Direct specification requires two arguments for efn.')
      return RNA.from_name_seq_db('cmdline', *args.cmd)


def process_command(*extra_args):
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', '--path', type=str)
  parser.add_argument('-kv', '--kekvault', type=str)
  parser.add_argument('cmd', type=str, nargs='*')

  parser.add_argument('--rnastructure-loc')
  parser.add_argument('--kekrna-loc')
  parser.add_argument('--viennarna-loc')
  parser.add_argument('--unafold-loc')
  parser.add_argument('--sparsemfefold-loc')

  parser.add_argument('-rh', '--rnastructure-harness', action='store_true')
  parser.add_argument('-rd', '--rnastructure-distribution', action='store_true')
  parser.add_argument('-m', '--rnark', action='store_true')
  parser.add_argument('-vd2', '--viennarna-d2', action='store_true')
  parser.add_argument('-vd3', '--viennarna-d3', action='store_true')
  parser.add_argument('-u', '--unafold', action='store_true')
  parser.add_argument('-k', '--kekrna', action='store_true')
  parser.add_argument('-smf', '--sparsemfefold', action='store_true')

  parser.add_argument('-f', '--fold', action='store_true')
  parser.add_argument('-e', '--energy', action='store_true')
  parser.add_argument('-s', '--subopt', type=str)
  parser.add_argument('-b', '--benchmark', type=str)

  args = parser.parse_args(sys.argv[1:] + list(*extra_args))

  if bool(args.fold) + bool(args.energy) + bool(args.subopt) + bool(args.benchmark) != 1:
    parser.error('Exactly one of --fold, --energy, --subopt, or --benchmark is required.')

  programs = []
  rnastructure_harness = RNAstructureHarness(args.kekrna_loc)
  if args.rnastructure_harness:
    programs.append(rnastructure_harness)
  if args.rnastructure_distribution:
    programs.append(RNAstructureDistribution(args.rnastructure_loc))
  if args.rnark:
    programs.append(Rnark(args.kekrna_loc))
  if args.viennarna_d2:
    programs.append(ViennaRNA(args.viennarna_loc, False))
  if args.viennarna_d3:
    programs.append(ViennaRNA(args.viennarna_loc, True))
  if args.unafold:
    programs.append(UNAFold(args.unafold_loc))
  if args.kekrna:
    programs.append(KekRNA(args.kekrna_loc))
  if args.sparsemfefold:
    programs.append(SparseMFEFold(args.sparsemfefold_loc))

  rna = parse_rna_from_args(parser, args)
  for program in programs:
    if args.fold:
      run_fold(program, rna)
    elif args.energy:
      run_efn(program, rna)
    elif args.subopt:
      run_suboptimal(program, rna, 6)  # TODO configurable
    elif args.benchmark:
      run_benchmark(program, args.benchmark, rnastructure_harness)

  for program in programs:
    program.close()


if __name__ == '__main__':
  process_command()
