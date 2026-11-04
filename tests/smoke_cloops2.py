#!/usr/bin/env python3
"""Opt-in cLoops2 smoke test on synthetic contacts only; never reads project samples."""
import argparse
import csv
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile


def synthetic_pairs(path):
    rng = random.Random(1739)
    contacts = []
    for anchor in (500000, 1500000, 3000000):
        for _ in range(80):
            contacts.append((anchor + rng.randrange(2000), anchor + 150000 + rng.randrange(2000)))
    for _ in range(600):
        left = rng.randrange(10000, 3500000)
        contacts.append((left, left + rng.randrange(25000, 1000000)))
    with path.open('w') as handle:
        handle.write('## pairs format v1.0.0\n#shape: upper triangle\n#chromsize: chr1 5000000\n')
        handle.write('#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n')
        for index, (left, right) in enumerate(sorted(contacts)):
            handle.write(f'toy{index}\tchr1\t{left}\tchr1\t{right}\t+\t-\tUU\n')


def run(command, root, env):
    return subprocess.run(command, cwd=root, env=env, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, timeout=120)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', default=shutil.which('cLoops2'))
    parser.add_argument('--check-parallel-failure', action='store_true')
    args = parser.parse_args()
    if not args.executable:
        parser.error('cLoops2 is not available; supply --executable /path/to/cLoops2')
    with tempfile.TemporaryDirectory(prefix='cloops2-synthetic-') as directory:
        root = Path(directory)
        env = os.environ.copy()
        env.update(MPLBACKEND='Agg', MPLCONFIGDIR=str(root / 'matplotlib'),
                   OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        pairs = root / 'synthetic.pairs'
        synthetic_pairs(pairs)
        pre = run([args.executable, 'pre', '-format', 'pairs', '-f', str(pairs),
                   '-o', str(root / 'pre'), '-c', 'chr1', '-cut', '20000', '-mcut', '2000000', '-p', '1'], root, env)
        if pre.returncode:
            raise RuntimeError(pre.stdout)
        base = [args.executable, 'callLoops', '-d', str(root / 'pre'), '-eps', '5000,10000',
                '-minPts', '10', '-hic', '-cut', '20000', '-mcut', '2000000', '-w', '-j']
        if args.check_parallel_failure:
            parallel = run(base + ['-p', '2', '-o', str(root / 'parallel')], root, env)
            print(f'Synthetic -p 2 exit status: {parallel.returncode}')
            if parallel.returncode and "'NoneType' object has no attribute 'info'" in parallel.stdout:
                print('Reproduced the macOS spawned-worker logger failure.')
            else:
                print('Parallel failure was not reproduced in this environment; see the original real-run log.')
                print(parallel.stdout[-2000:])
        serial = run(base + ['-p', '1', '-o', str(root / 'serial')], root, env)
        if serial.returncode:
            raise RuntimeError(serial.stdout)
        with (root / 'serial_loops.txt').open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not reader.fieldnames or reader.fieldnames[0] != 'loopId':
                raise RuntimeError('Native cLoops2 output has an invalid header')
            count = sum(1 for _ in reader)
        print(f'Serial cLoops2 CLI smoke test passed on 840 synthetic contacts; output rows={count}.')
        print('No biological recovery or real-data analysis is validated by this synthetic test.')


if __name__ == '__main__':
    main()
