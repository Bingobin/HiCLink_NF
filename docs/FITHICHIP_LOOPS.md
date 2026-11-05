# FitHiChIP branch of the reusable loop subworkflow

The branch has four linear steps: prepare raw matrices, merge peaks, call loops, and run the official differential script. It is an optional caller within `HICHIP_LOOPS`, not a replacement for `main.nf`.

## Install the Docker runtime

```bash
bash bin/install_fithichip_docker.sh
```

The official `aylab/fithichip:latest` image retrieved on 2026-11-05 contains older scripts, lacks Python 3/cooler and lacks R data.table. The supplied Dockerfile **extends the official image** pinned to digest `sha256:bbb2cbcfc4bc61123d5100b111b2a648100da5fe8d24236ca0e4a3b7d8162fb4`, installs the missing dependencies, and copies the unmodified official source at commit `0ea1ac21be870908c672316ffbb630189dc6fae2`. No caller/statistical algorithm is patched. The result is the local image `hiclink/fithichip:11-0ea1ac2`.

The runtime is `linux/amd64`; Apple Silicon runs it through Docker's compatibility layer. It does not install anything into the host R library. Docker must be running. The build directory is retained and printed by the installer.

The requested shell interface is retained:

```bash
bash bin/FitHiChIP_Docker.sh -C /absolute/path/sample.config.txt
```

This is a small workflow adapter, **not a verbatim copy of the upstream Docker wrapper**. The upstream wrapper checks host installations of HiC-Pro/R/MACS2 before launching Docker. This adapter instead runs dependencies inside Docker, mounts named inputs read-only and a new output directory read-write, and invokes the official `FitHiChIP_HiCPro.sh`. `--image`, `--platform`, `--cpus` and `--memory-bytes` can override its runtime settings. `COOL`, `PeakFile`, `ChrSizeFile`, `OutDir` and the scientific parameters follow the upstream config format. `OutDir` must be new. Do not also enable Nextflow's generic Docker profile for these wrapper tasks.

## Run from published inputs

```bash
cp conf/hichip_loops_fithichip.example.config conf/hichip_loops_fithichip.local.config
# Edit the paths and existing host cooler/Python environment.
nextflow run hichip_loops.nf \
  -c conf/hichip_loops_fithichip.local.config -profile hichip_loops_local
```

Use a new output version such as `v2_fithichip_10kb`; do not reuse the historical cLoops2/Mustache output directory. Unchanged retries can use `-resume`. The shared helper/settings changed, so the original signature guard intentionally rejects an old output directory rather than mixing versions.

Each sample still needs the existing common `tuple(meta, pairs, matrix)` input. FitHiChIP does not read the pairs records. It reads raw matrix counts, not ICE-weighted counts. Host preparation needs cooler and its Python dependencies. Native FitHiChIP/R dependencies are provided by Docker.

The standalone entry discovers exact `06_hichip/peaks/ID.peaks.bed` files, or accepts the explicit `peak_file` CSV column. Missing treatment peaks fail; control peaks are not silently substituted. Supplied paths in `peak_file` should be absolute. In the parent-workflow interface, set `meta.peak_file` to a staged/accessible BED path. Peaks must use the same genome and chromosome naming as the matrix.

## Steps and outputs

| Step | Input | Output / meaning |
| --- | --- | --- |
| `FITHICHIP_PREPARE` | Per-sample cool/mcool + peaks | `fithichip/inputs/ID.fit_inputs/`: raw-count `contacts.cool`, chromosome sizes and normalized BED peaks; configured chromosomes only; exact/finer integer-multiple resolution allowed |
| `FITHICHIP_UNION_PEAKS` | Control and treatment peaks | One merged union BED per comparison; both samples use exactly this peak universe |
| `FITHICHIP_LOOP_CALL` | Each raw matrix + comparison's union peaks | Independent FitHiChIP output for each role; all native files retained, plus an unthresholded `ID.interactions_FitHiC.bed` at the result root |
| `FITHICHIP_LOOP_DIFF` | **Unthresholded** control and treatment FitHiChIP tables | Official `Imp_Scripts/DiffAnalysisHiChIP.r` output plus an annotated all-tested table and `comparison.json` |

A control shared by multiple comparisons is refitted with each comparison's union peaks; those results must not be interchanged. Samples without a comparison, such as HMGA2-OE-HA, are called individually using their own peaks. There is no cross-IP differential comparison. Original matrices and peaks are not modified or overwritten. A finer source is coarsened into a new temporary file; weights are not copied or used in the prepared raw matrix.

The existing core signature is retained: `HICHIP_LOOPS(samples, comparisons, cfg)`. Additional emitted channels are `fithichip_inputs`, `fithichip_calls`, and `fithichip_diff`. Call metadata adds `run_id`, `role`, and, for paired calls, `comparison`.

## Parameters

| Parameter | Default |
| --- | --- |
| `hichip_loop_callers` | Set to `fithichip`, or combine with existing callers |
| `fithichip_loop_image` / `fithichip_loop_platform` | `hiclink/fithichip:11-0ea1ac2` / `linux/amd64` |
| `fithichip_loop_cpus` / `fithichip_loop_memory_gb` | 2 / 6; Docker resource limits, not measured requirements |
| `fithichip_loop_interaction_type` | 3, peak-to-all |
| `fithichip_loop_bias_type` | 1, coverage-bias regression |
| `fithichip_loop_p2p_background` | 0, loose peak-to-all background |
| `fithichip_loop_qvalue` | 0.01, single-sample FitHiChIP threshold |
| `fithichip_loop_background_fdr` | 1.0, differential-model background universe |
| `fithichip_loop_diff_fdr` | 0.05, edgeR BH FDR threshold |
| `fithichip_loop_fold_change` | 2, absolute fold change; upstream uses log2(2)=1 |
| `fithichip_loop_bcv` | 0.4, assumed BCV for the current single-library comparisons |

Resolution, chromosome selection and distance bounds use the existing loop settings. The current module interface specifies **one control and one treatment library per comparison**. The official R script is invoked with `Control:Treatment` and `1:1`; positive logFC is treatment-enriched. Do not claim replicate-aware inference or treat separate cell lines as replicates. Multi-replicate condition grouping is not implemented by this sample-pair interface.

`differential.annotated.tsv` preserves all native edgeR columns and adds `significant` and `direction`. Significant requires differential FDR/fold-change thresholds **and** FitHiChIP significance in at least one sample, matching the official filtered output. Non-significant rows are explicitly labelled. Native edgeR and single-sample FitHiChIP FDRs remain separate.

Matching independent ChIP BAM/bedGraph inputs are not supplied by this module, so official 1D-occupancy HD/ND/LD stratification is **not performed**. Do not claim that differential HiChIP signal is solely a 3D change. No RPGC BigWig is substituted for ChIP counts. Union peaks avoid different anchor definitions, but do not remove all occupancy/IP/depth effects.

## Checks

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
nextflow lint hichip_loops.nf subworkflows/local/hichip_loops.nf modules/local/fithichip_loops.nf
python tests/make_fithichip_fixture.py /tmp/fithichip-synthetic
```

The fixture creates synthetic contacts and peaks only. It can be used with the standalone entry and an isolated output/work directory to test real Docker execution; it is not a biological benchmark. Production HMGA2 calling is not automatically started by installation or tests.

References: [official Docker wrapper](https://github.com/ay-lab/FitHiChIP/blob/master/FitHiChIP_Docker.sh), [official differential analysis](https://ay-lab.github.io/FitHiChIP/html/usage/DiffLoops.html), and [official source](https://github.com/ay-lab/FitHiChIP).
