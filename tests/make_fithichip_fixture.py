#!/usr/bin/env python3
"""Generate reproducible synthetic contacts for an optional end-to-end Docker smoke test."""
from pathlib import Path
import argparse
import cooler
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=Path)
args = parser.parse_args()
root = args.output.resolve()
for directory in ['03_pairs', '04_matrices', '06_hichip/peaks', 'pipeline_info']:
    (root / directory).mkdir(parents=True, exist_ok=True)
n = 350
bins = pd.DataFrame({'chrom':['chr1']*n, 'start':np.arange(n)*10000, 'end':(np.arange(n)+1)*10000})
a, b = np.triu_indices(n)
distance=b-a
base = np.maximum(2, np.rint(30/np.sqrt(distance+1))).astype(int)
for sample, multiplier in [('CTRL',1),('KD',300),('EXTRA',1)]:
    rng=np.random.default_rng(17)
    counts=rng.poisson(base)+1
    spike=(a % 23 == 0)&(distance==15)
    counts[spike] += 30*multiplier
    pixels=pd.DataFrame({'bin1_id':a,'bin2_id':b,'count':counts})
    cooler.create_cooler(str(root/f'04_matrices/{sample}.10000.cool'),bins,pixels)
    (root/f'03_pairs/{sample}.valid.pairs').write_text('# Synthetic placeholder: the FitHiChIP-only branch never reads PET records.\n')
    peaks=''.join(f'chr1\t{i*10000}\t{i*10000+5000}\n' for i in range(0,n,3))
    (root/f'06_hichip/peaks/{sample}.peaks.bed').write_text(peaks)
(root/'pipeline_info/samplesheet.validated.csv').write_text('ID,assay,control\nCTRL,hichip,\nKD,hichip,CTRL\nEXTRA,hichip,\n')
print(root)
