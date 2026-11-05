"""FitHiChIP adapters: synthetic inputs only, no biological data or caller execution."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import cooler
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import fithichip as FIT
import hichip_loops as LOOPS
from test_hichip_loops import settings


class FitHiChIPTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_union_merges_overlap_and_preserves_sources(self):
        a, b, out = [self.root / name for name in ['a.bed', 'b.bed', 'union.bed']]
        a.write_text('chr1\t0\t100\nchr2\t0\t10\n')
        b.write_text('chr1\t90\t150\nchr1\t300\t400\n')
        FIT.merge_peaks([a, b], out, ['chr1'])
        self.assertEqual(out.read_text(), 'chr1\t0\t150\nchr1\t300\t400\n')
        self.assertIn('chr2', a.read_text())

    def test_prepare_keeps_raw_counts_and_selected_chromosomes(self):
        matrix = self.root / 'source.cool'
        bins = pd.DataFrame({'chrom':['chr1','chr1','chr2'], 'start':[0,10000,0], 'end':[10000,20000,10000]})
        pixels = pd.DataFrame({'bin1_id':[0,0,0,2], 'bin2_id':[0,1,2,2], 'count':[4,7,3,9]})
        cooler.create_cooler(str(matrix), bins, pixels)
        peaks = self.root / 'peaks.bed'
        peaks.write_text('chr1\t100\t200\nchr2\t100\t200\n')
        FIT.prepare(matrix, peaks, self.root / 'prepared', settings(), 1)
        clr = cooler.Cooler(str(self.root / 'prepared/contacts.cool'))
        self.assertEqual(clr.chromnames, ['chr1'])
        self.assertEqual(clr.info['sum'], 11)
        self.assertEqual(clr.matrix(balance=False)[:].tolist(), [[4,7],[7,0]])
        self.assertNotIn('weight', clr.bins().columns)
        self.assertEqual(cooler.Cooler(str(matrix)).info['sum'], 23)

    def test_fit_only_discovery_requires_exact_peaks_and_tracks_them(self):
        for directory in ['03_pairs','04_matrices','06_hichip/peaks']:
            (self.root / directory).mkdir(parents=True)
        for name in ['03_pairs/S1.valid.pairs.gz','04_matrices/S1.mcool']:
            (self.root / name).touch()
        sheet = self.root / 'samples.csv'
        sheet.write_text('ID,assay\nS1,hichip\n')
        cfg = settings()
        cfg.update(callers=['fithichip'], outdir=str(self.root/'published'),
                   fithichip={'interaction_type':3,'fold_change':2,'bcv':0.4})
        with self.assertRaisesRegex(ValueError, 'Missing peaks for S1'):
            LOOPS.discover(sheet,self.root,cfg,self.root/'missing')
        peak = self.root/'06_hichip/peaks/S1.peaks.bed'
        peak.write_text('chr1\t0\t100\n')
        LOOPS.discover(sheet,self.root,cfg,self.root/'manifest')
        signature=json.loads((self.root/'manifest/loop_inputs_signature.json').read_text())
        self.assertIn(str(peak.resolve()),[r['path'] for r in signature['input_stats']])

    def test_configuration_uses_coverage_bias_and_original_raw_matrix(self):
        cfg = settings()
        cfg['fithichip'] = {'interaction_type':3,'p2p_background':0,'bias_type':1,'qvalue':0.01}
        FIT.write_config(self.root,self.root/'union.bed',self.root/'new-result','S1',cfg,self.root/'config.txt')
        contents=(self.root/'config.txt').read_text()
        for value in ['BiasType=1','UseP2PBackgrnd=0','OverWrite=0','PREFIX=S1']:
            self.assertIn(value, contents)


if __name__ == '__main__':
    unittest.main()
