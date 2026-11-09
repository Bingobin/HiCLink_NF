"""Regression for ordered Nextflow config evaluation; no analysis task is started."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')
if not NEXTFLOW and (Path.home() / 'bin/nextflow').is_file():
    NEXTFLOW = str(Path.home() / 'bin/nextflow')


def quote(value):
    return "'" + str(value).replace('\\', '\\\\').replace("'", "\\'") + "'"


@unittest.skipUnless(NEXTFLOW, 'Nextflow is required for configuration-only regression checks')
class LoopRuntimeConfigTests(unittest.TestCase):
    def test_local_concurrency_is_resource_driven_not_a_fixed_task_count(self):
        with tempfile.TemporaryDirectory(prefix='loop-scheduling-') as directory:
            root = Path(directory)
            driver = root / 'job.config'
            driver.write_text(
                f"includeConfig {quote(ROOT / 'conf/hichip_loops.config')}\n"
                f"includeConfig {quote(ROOT / 'conf/hichip_loops_runtime.config')}\n"
            )
            result = subprocess.run(
                [NEXTFLOW, '-c', str(driver), 'config', str(ROOT), '-profile', 'hichip_loops_local', '-flat'],
                cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env={**os.environ, 'NXF_OFFLINE': 'true'}, timeout=60, check=True,
            )
            self.assertIn("process.executor = 'local'", result.stdout)
            self.assertIn("process.'withLabel:hichip_loop_analysis'.memory = {", result.stdout)
            self.assertIn("process.'withName:.*CLOOPS_CALL'.cpus = {", result.stdout)
            for line in result.stdout.splitlines():
                key = line.partition(' = ')[0]
                self.assertNotIn('maxForks', key)
                self.assertNotIn('queueSize', key)
                self.assertNotIn(key, ('executor.cpus', 'executor.memory'))

    def test_ordered_config_rebinds_all_report_paths_after_job_parameters(self):
        with tempfile.TemporaryDirectory(prefix='loop-config-') as directory:
            root = Path(directory)
            driver = root / 'job.config'
            output = root / 'actual_analysis'
            driver.write_text(
                f"includeConfig {quote(ROOT / 'conf/hichip_loops.config')}\n"
                'params {\n'
                f"  outdir = {quote(root / 'outer_results')}\n"
                f"  hichip_loop_outdir = {quote(output)}\n"
                '  hichip_loop_cpus = 6\n'
                '}\n'
                f"includeConfig {quote(ROOT / 'conf/hichip_loops_runtime.config')}\n"
            )
            result = subprocess.run(
                [NEXTFLOW, '-c', str(driver), 'config', str(ROOT), '-profile', 'hichip_loops_local', '-flat'],
                cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env={**os.environ, 'NXF_OFFLINE': 'true'}, timeout=60, check=True,
            )
            for scope, name in [('report', 'execution_report.html'), ('timeline', 'execution_timeline.html'),
                                ('trace', 'execution_trace.txt'), ('dag', 'pipeline_dag.html')]:
                self.assertIn(f"{scope}.file = '{output}/pipeline_info/{name}'", result.stdout)
            self.assertIn('params.hichip_loop_cpus = 6', result.stdout)
            self.assertIn("process.'withName:.*CLOOPS_CALL'.cpus", result.stdout)
            self.assertIn("process.'withLabel:hichip_loop_analysis'.cpus = {", result.stdout)
            self.assertNotIn("file = 'null/", result.stdout)
            self.assertFalse(output.exists(), 'Configuration inspection must not create analysis outputs')


if __name__ == '__main__':
    unittest.main()
