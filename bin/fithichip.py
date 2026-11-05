#!/usr/bin/env python3
"""Small input/command adapter for the official FitHiChIP Docker workflow."""
import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess

DEFAULT_IMAGE = 'hiclink/fithichip:11-0ea1ac2'


def merge_peaks(sources, output, chroms=None):
    """Merge BED intervals without changing the original peak files."""
    intervals = []
    for source in sources:
        opener = gzip.open if str(source).endswith('.gz') else open
        with opener(source, 'rt') as handle:
            for line in handle:
                if line.startswith(('#', 'track', 'browser')) or not line.strip():
                    continue
                chrom, start, end = line.split()[:3]
                if chroms is None or chrom in chroms:
                    intervals.append((chrom, int(start), int(end)))
    merged = []
    for chrom, start, end in sorted(intervals):
        if merged and chrom == merged[-1][0] and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    if not merged:
        raise ValueError('No peaks overlap the configured chromosomes')
    with Path(output).open('w') as handle:
        for interval in merged:
            handle.write('\t'.join(map(str, interval)) + '\n')


def prepare(matrix, peaks, output, cfg, cpus):
    import cooler
    import numpy as np
    from hichip_loops import matrix_source

    output.mkdir()
    uri, clr = matrix_source(matrix, cfg['resolution'])
    if clr.binsize != cfg['resolution']:
        coarse = output / 'coarse.cool'
        cooler.coarsen_cooler(uri, str(coarse), cfg['resolution'] // clr.binsize,
                             chunksize=1_000_000, nproc=cpus)
        clr = cooler.Cooler(str(coarse))
    bins = clr.bins()[:][['chrom', 'start', 'end']]
    keep = bins.chrom.isin(cfg['chroms']).to_numpy()
    index = np.full(len(bins), -1, dtype=np.int64)
    index[keep] = np.arange(keep.sum())

    def pixels():
        for offset in range(0, clr.info['nnz'], 1_000_000):
            frame = clr.pixels()[offset:offset + 1_000_000].copy()
            frame = frame[keep[frame.bin1_id] & keep[frame.bin2_id]].copy()
            frame['bin1_id'] = index[frame.bin1_id]
            frame['bin2_id'] = index[frame.bin2_id]
            yield frame

    cooler.create_cooler(str(output / 'contacts.cool'), bins[keep].reset_index(drop=True),
                         pixels(), ordered=True, dtypes={'count': np.int64})
    with (output / 'chrom.sizes').open('w') as handle:
        for chrom in cfg['chroms']:
            handle.write(f'{chrom}\t{int(clr.chromsizes[chrom])}\n')
    merge_peaks([peaks], output / 'peaks.bed', cfg['chroms'])
    print(f'Prepared raw {cfg["resolution"]}-bp contacts and peaks in {output}')


def write_config(inputs, peaks, output_dir, sample, cfg, target):
    fit = cfg['fithichip']
    values = dict(COOL=str((inputs / 'contacts.cool').resolve()),
                  PeakFile=str(peaks.resolve()), ChrSizeFile=str((inputs / 'chrom.sizes').resolve()),
                  OutDir=str(output_dir.resolve()), PREFIX=sample, CircularGenome=0,
                  IntType=fit['interaction_type'], BINSIZE=cfg['resolution'],
                  LowDistThr=cfg['min_distance'], UppDistThr=cfg['max_distance'],
                  UseP2PBackgrnd=fit['p2p_background'], BiasType=fit['bias_type'],
                  MergeInt=1, QVALUE=fit['qvalue'], OverWrite=0)
    target.write_text(''.join(f'{key}={value}\n' for key, value in values.items()))


def docker_run(args, command, inputs, output):
    """Mount only named inputs read-only and the new result directory read-write."""
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cmd = ['docker', 'run', '--rm', '--platform', args.platform, '--network', 'none',
           '--cpus', str(args.cpus), '--memory', str(args.memory_bytes),
           '--env', 'OMP_NUM_THREADS=1', '--env', 'OPENBLAS_NUM_THREADS=1',
           '--env', 'MKL_NUM_THREADS=1', '--workdir', '/FitHiChIP',
           '--mount', f'type=bind,src={output},dst=/output']
    for source, target in inputs:
        cmd += ['--mount', f'type=bind,src={Path(source).resolve()},dst={target},readonly']
    cmd += ['--entrypoint', command[0], args.image, *command[1:]]
    (output / 'docker.command.txt').write_text(shlex.join(cmd) + '\n')
    image_id = subprocess.check_output(['docker', 'image', 'inspect', args.image, '--format', '{{.Id}}'], text=True)
    (output / 'docker.image.txt').write_text(image_id)
    print(shlex.join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def run_call(args):
    values = dict(line.strip().split('=', 1) for line in args.config.read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#'))
    output = Path(values['OutDir'])
    if output.exists():
        raise ValueError(f'Refusing to reuse an existing result directory: {output}')
    output.mkdir(parents=True)
    chrom_sizes = Path(values['ChrSizeFile'])
    inputs = []
    for key in ['COOL', 'PeakFile', 'ChrSizeFile']:
        source = Path(values[key]).resolve()
        if not source.is_file():
            raise FileNotFoundError(f'{key} must name an existing input file: {source}')
        target = '/input/' + key + ''.join(source.suffixes)
        inputs.append((source, target))
        values[key] = target
    values['OutDir'] = '/output/native'
    container_config = output / 'container.config.txt'
    container_config.write_text(''.join(f'{k}={v}\n' for k, v in values.items()))
    docker_run(args, ['/bin/bash', '/FitHiChIP/FitHiChIP_HiCPro.sh',
                      '-C', '/output/container.config.txt'], inputs, output)
    tables = list((output / 'native').rglob(values['PREFIX'] + '.interactions_FitHiC.bed'))
    if len(tables) != 1:
        raise ValueError(f'Expected one unthresholded FitHiChIP table, found {len(tables)}; inspect the task log')
    shutil.copy2(tables[0], output / (values['PREFIX'] + '.interactions_FitHiC.bed'))
    shutil.copy2(args.config, output / 'host.config.txt')
    shutil.copy2(chrom_sizes, output / 'chrom.sizes')


def run_diff(args, cfg):
    fit = cfg['fithichip']
    command = ['Rscript', '/FitHiChIP/Imp_Scripts/DiffAnalysisHiChIP.r',
               '--AllLoopList', '/input/control.bed:/input/treatment.bed',
               '--ChrSizeFile', '/input/chrom.sizes', '--OutDir', '/output/native',
               '--CategoryList', 'Control:Treatment', '--ReplicaCount', '1:1',
               '--ReplicaLabels1', 'R1', '--ReplicaLabels2', 'R1',
               '--FDRThr', str(fit['qvalue']), '--BackgroundFDRThr', str(fit['background_fdr']),
               '--FoldChangeThr', str(fit['fold_change']), '--DiffFDRThr', str(fit['diff_fdr']),
               '--bcv', str(fit['bcv'])]
    docker_run(args, command, [(args.control, '/input/control.bed'),
                               (args.treatment, '/input/treatment.bed'),
                               (args.chrom_sizes, '/input/chrom.sizes')], args.output)
    native = args.output / 'native/EdgeR_Loops_ALL/Loops_EdgeR_Default.bed'
    if not native.is_file():
        raise ValueError('Missing official edgeR output; inspect the task log')
    with native.open() as handle, (args.output / 'differential.annotated.tsv').open('w') as out:
        reader = csv.DictReader(handle, delimiter='\t')
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames + ['significant', 'direction'], delimiter='\t')
        writer.writeheader()
        for row in reader:
            row['significant'] = int(float(row['FDR']) <= fit['diff_fdr'] and
                                     abs(float(row['logFC'])) >= math.log2(fit['fold_change']) and
                                     (int(row['Control_SigRepl']) > 0 or int(row['Treatment_SigRepl']) > 0))
            row['direction'] = ('treatment_enriched' if float(row['logFC']) > 0 else 'control_enriched') if row['significant'] else 'not_significant'
            writer.writerow(row)
    (args.output / 'comparison.json').write_text(json.dumps(args.comparison, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--matrix', type=Path, required=True)
    p.add_argument('--peaks', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--cpus', type=int, default=1)
    p.add_argument('--settings', type=json.loads, required=True)
    p = sub.add_parser('union')
    p.add_argument('sources', type=Path, nargs='+')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('config')
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--peaks', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--settings', type=json.loads, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run-call')
    p.add_argument('-C', '--config', type=Path, required=True)
    p = sub.add_parser('diff')
    p.add_argument('--control', type=Path, required=True)
    p.add_argument('--treatment', type=Path, required=True)
    p.add_argument('--chrom-sizes', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--comparison', type=json.loads, required=True)
    p.add_argument('--settings', type=json.loads, required=True)
    for name in ['run-call', 'diff']:
        p = sub.choices[name]
        p.add_argument('--image', default=os.environ.get('FITHICHIP_IMAGE', DEFAULT_IMAGE))
        p.add_argument('--platform', default='linux/amd64')
        p.add_argument('--cpus', type=int, default=2)
        p.add_argument('--memory-bytes', type=int, default=6 * 1024 ** 3)
    a = parser.parse_args()
    if a.command == 'prepare':
        prepare(a.matrix, a.peaks, a.output, a.settings, a.cpus)
    elif a.command == 'union':
        merge_peaks(a.sources, a.output)
    elif a.command == 'config':
        write_config(a.inputs, a.peaks, a.output_dir, a.sample, a.settings, a.output)
    elif a.command == 'run-call':
        run_call(a)
    elif a.command == 'diff':
        run_diff(a, a.settings)


if __name__ == '__main__':
    main()
