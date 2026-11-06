# HiCLink_NF

`HiCLink_NF` is a Nextflow DSL2 pipeline for unified Hi-C and HiChIP analysis.

The original `hic_cool.nf` and `hic_proc.nf` scripts are retained as legacy entry points. The new workflow starts from `main.nf`.

## Workflow

```text
FASTQ
  ├─ FastQC and fastp independently for each sequencing run
  ├─ merge cleaned R1/R2 FASTQs by library ID
  ├─ BWA-MEM alignment (-SP)
  ├─ pairtools parse, sort, deduplication, filtering and QC
  ├─ balanced single- and multi-resolution cooler matrices
  ├─ shared: contact coverage and matrix sparsity QC
  ├─ Hi-C: expected-cis, compartments, insulation/TADs and cooltools dots
  └─ HiChIP: RPGC BigWig, supplied or MACS2 peaks and HOMER annotation
```

A single samplesheet can contain both `hic` and `hichip` libraries. The pipeline retains `condition`, `replicate` and `control` metadata. An optional cLoops2/Mustache/FitHiChIP loop subworkflow and published-file entry point are available separately; replicate-aware inference is not implemented by that sample-pair module.

HiChIP matrices contain target-enrichment bias. Compartment and TAD analysis therefore runs on Hi-C samples by default. Use `--global_features_all_assays true` only when there is a specific methodological reason to calculate these features for HiChIP, and interpret them separately from unbiased Hi-C maps.

## Samplesheet

Use `assets/samplesheet.example.csv` as the template.

| Column | Required | Description |
|---|---:|---|
| `ID` | Yes | Physical library ID containing letters, numbers, `.`, `_` or `-`; repeat it for additional sequencing runs of that same library |
| `R1`, `R2` | Yes | Paired-end FASTQ files; relative paths are resolved from the samplesheet directory |
| `assay` | Yes | `hic` or `hichip` |
| `condition` | No | Biological condition; defaults to the sample ID |
| `replicate` | No | Biological replicate number; defaults to 1 |
| `control` | No | Sample ID of a matched control; used by the optional HiChIP loop entry point |
| `peak_file` | HiChIP only | External ChIP-seq or consensus peaks; MACS2 is used when this field is empty |

Example:

```csv
ID,R1,R2,assay,condition,replicate,control,peak_file
HIC_CTRL_R1,/data/hic_ctrl_R1.fastq.gz,/data/hic_ctrl_R2.fastq.gz,hic,CTRL,1,,
HICHIP_CTRL_R1,/data/hichip_ctrl_R1.fastq.gz,/data/hichip_ctrl_R2.fastq.gz,hichip,CTRL,1,,/data/H3K27ac.peaks.bed
HICHIP_TREAT_R1,/data/hichip_R1.fastq.gz,/data/hichip_R2.fastq.gz,hichip,TREAT,1,HICHIP_CTRL_R1,/data/H3K27ac.peaks.bed
HICHIP_TREAT_R1,/data/hichip_run2_R1.fastq.gz,/data/hichip_run2_R2.fastq.gz,hichip,TREAT,1,HICHIP_CTRL_R1,/data/H3K27ac.peaks.bed
```

### Multiple sequencing runs of one library

Use one row per paired FASTQ input and repeat the same `ID`. All rows of an ID must agree on `assay`, `condition`, `replicate`, `control` and normalized `peak_file`. Reusing a FASTQ path within a library is rejected to avoid double counting. Different biological replicates and independently prepared libraries must retain different IDs, even if they share a condition.

The validator adds an output-only `run_id` column. A multi-run library receives `ID.run001`, `ID.run002`, etc., in samplesheet order; a single-row library retains `ID`. No extra input column or merge parameter is required. FastQC and fastp run independently per row, using these distinct names in reports and MultiQC. After all runs finish cleaning, `MERGE_FASTQ` concatenates their gzip files in matching R1/R2 order, then runs BWA, pairtools and matrix generation once per library. Single-run inputs use symlinks instead of copying their FASTQs. The original raw/clean input files are not modified.

`pipeline_info/samplesheet.validated.csv` maps raw input paths to generated run IDs, and `pipeline_info/sequencing_runs/ID.sequencing_runs.tsv` records the clean FASTQs used for each merge. Merged FASTQs remain in the Nextflow work directory; `save_trimmed` still publishes the per-run fastp outputs. Keep original sequencer read identifiers so reads from different runs remain distinguishable.

Final BAM, pairs, cool/mcool and peak outputs retain the original library `ID`. Pairtools duplicate handling takes place after merging, so `skip_pairtools_dedup` applies to the combined library. Feature restarts and the standalone loop entry point collapse repeated run rows to one library and reject conflicting metadata, without accessing FASTQs. Additional sequencing does not create additional biological replicates.

## Reference files

Required reference files:

- A reference FASTA with classic BWA index files beside it: `.amb`, `.ann`, `.bwt`, `.pac` and `.sa`.
- A two-column chromosome sizes file without a header.
- An optional four-column, header-containing phasing track for compartment analysis. Its resolution must match `compartment_resolution`, and `phasing_track_column` identifies its signal column, for example `GC`.
- An effective genome size is required when `save_peak_bigwig` is enabled. `genome: hg38` or `GRCh38` automatically uses `2913022398`; `effective_genome_size` can override it and is required for custom genomes.

Reference preparation example:

```bash
bwa index /reference/hg38.fa
samtools faidx /reference/hg38.fa
cut -f1,2 /reference/hg38.fa.fai > /reference/hg38.chrom.sizes
```

## Installation

Nextflow `>=24.10.0` and Java are required. Core bioinformatics dependencies can be installed automatically with Conda or Mamba using `envs/hiclink.yml`.

FitHiChIP is provided only by the optional HiChIP loop subworkflow, using its Docker launcher and `fithichip_loop_*` settings. The main workflow does not run FitHiChIP. See [HiChIP loop documentation](docs/HICHIP_LOOPS.md); the legacy main-workflow FitHiChIP process, parameters and container profiles have been removed.

## Usage

Copy and edit the example configuration and samplesheet:

```bash
cp conf/params.example.yml conf/params.local.yml
cp assets/samplesheet.example.csv samplesheet.csv
```

Run with tools already available in `PATH`:

```bash
nextflow run main.nf \
  -profile standard \
  -params-file conf/params.local.yml \
  --input samplesheet.csv
```

Run with the core Conda environment:

```bash
nextflow run main.nf \
  -profile conda \
  -params-file conf/params.local.yml \
  --input samplesheet.csv
```

Run the main workflow on SLURM:

```bash
nextflow run main.nf \
  -profile conda,slurm \
  -params-file conf/params.local.yml \
  --input samplesheet.csv
```

For a HiChIP sample without `peak_file`, set `--macs2_genome_size` to `hs`, `mm` or an effective genome size. The workflow converts primary HiChIP alignments to BED before MACS2 peak calling because proximity-ligation pairs are not ordinary short-insert `BAMPE` fragments. MACS2 runs with `--nomodel`; `macs2_extsize` controls the fixed extension length and defaults to 147 bp.

### Disable loop calling

Append `--run_loops false`, or set `run_loops: false` in the parameter YAML. The default is `true`. The main workflow continues per-run QC, library merging, alignment, pairs, matrices and other enabled features, but skips cooltools dots. Peak calling and annotation remain independently controlled.

The same switch disables all callers and differential-loop branches of the optional HiChIP loop workflow. Its standalone entry point exits without input discovery or analysis when disabled. With `run_loops: true`, the caller-specific settings (`run_hic_loops` in the main workflow and `hichip_loop_callers` in the loop workflow) still apply. HiChIP FitHiChIP calling is never duplicated in the main workflow.

### Restart feature analysis from published results

If `work/` was deleted, start a new feature run using the published results rather than relying on `-resume`:

```bash
nextflow run main.nf \
  -profile standard \
  -params-file conf/params.local.yml \
  --start_from features \
  --previous_results /analysis/previous_run/results \
  --input /analysis/previous_run/results/pipeline_info/samplesheet.validated.csv \
  --outdir results_features
```

`start_from` defaults to `fastq`, which runs the full workflow. In `features` mode, the workflow validates sample metadata and reads `04_matrices/ID.mcool` from `previous_results`. It skips FastQC, fastp, alignment, pairtools, matrix generation, BigWig and MultiQC. FASTQ files, the reference FASTA and BWA indexes are not accessed. The samplesheet needs only `ID` and `assay`; the existing full samplesheet can also be used. Keep sample IDs identical to the published filenames.

For HiChIP with `run_hichip_peaks: true`, peak inputs are selected in this order:

1. `previous_results/06_hichip/peaks/ID.peaks.bed` (reuse without MACS2).
2. The samplesheet's `peak_file` (including `.narrowPeak` or `.gz`; no BAM required).
3. `previous_results/02_preprocessing/alignment_bam/ID.peaks.bam` and `.bam.bai` (run MACS2; requires that `save_alignment_bam` was enabled in the original run).

Missing matrices or required peak inputs produce an error naming the sample. For matrix features only, set `--run_hichip_peaks false --run_peak_annotation false`. When running peak annotation, enable `run_hichip_peaks` even if peaks are already available: this step imports the peak file instead of calling MACS2. FitHiChIP uses the separate loop entry point. Structural HiChIP features still require `--global_features_all_assays true`; the default Hi-C/HiChIP routing is unchanged.

Keep `genome`, `chrom_sizes`, resolutions and the optional compartment phasing track consistent with the original matrices. The requested resolutions and balancing weights must already exist in each `.mcool`; feature mode does not rebuild or rebalance it. Published files must be real copies or working links; links into a deleted work directory cannot be recovered by this mode. A separate `outdir` preserves the old report and feature outputs. Subsequent retries of the new feature run can use `-resume` as long as its new work directory is retained.

If one sample's matrix build was interrupted but its valid pairs were published, rebuild that sample independently before starting features:

```bash
conda activate hiclink-nf
bash bin/rebuild_cooler.sh A549.RAD21.SCR \
  results/03_pairs/A549.RAD21.SCR.valid.pairs.gz \
  /reference/hg38.chrom.sizes \
  results/04_matrices
```

The standalone Bash script runs `cload pairs`, `balance` and `zoomify --balance`, matching the workflow defaults (hg38, 16 CPUs, 5 kb base resolution, the standard eight resolutions, and `ignore_diags=2`). If the original run used different settings, set `HICLINK_GENOME`, `HICLINK_CPUS`, `HICLINK_BASE_RESOLUTION`, `HICLINK_RESOLUTIONS` and/or `HICLINK_IGNORE_DIAGS` before invoking it. Both `.cool` and `.mcool` are published only after the commands succeed. Existing output files are never overwritten; use a separate output directory if necessary. Partial builds are retained in the printed temporary directory on failure.

### Optional cLoops2/Mustache/FitHiChIP loops

For the Docker-based FitHiChIP branch, install with `bash bin/install_fithichip_docker.sh` and use `conf/hichip_loops_fithichip.example.config`. It prepares raw counts, shares the control/treatment union peaks, calls loops and runs the official differential R script. See [FitHiChIP setup and outputs](docs/FITHICHIP_LOOPS.md). The `bash bin/FitHiChIP_Docker.sh -C CONFIG` interface is also available for standalone calling. The original FASTQ entry point is unchanged.

A reusable `HICHIP_LOOPS` subworkflow is available in `subworkflows/local/hichip_loops.nf`. It is not yet called by `main.nf`. Its separate `hichip_loops.nf` entry point reads an existing `samplesheet.validated.csv`, finds matching published pairs/cool/mcool by exact sample ID, and derives directed comparisons from `control`. FASTQ and reference files are not accessed.

```bash
cp conf/hichip_loops_job.example.config conf/hichip_loops_job.local.config
# Edit input, previous_results, output paths, and the optional existing environment.
nextflow run hichip_loops.nf -c conf/hichip_loops_job.local.config \
  -profile hichip_loops_local
```

See [the subworkflow documentation](docs/HICHIP_LOOPS.md) for the channel contract, future integration example, defaults, output directions, caching, and validation limits. This entry point does not install loop dependencies. With one library per group, differential results are exploratory, not replicate-aware inference. No real-data loop analysis was run during implementation.

The ordered job config loads base settings, job overrides, then runtime report bindings. This keeps reports beside loop results instead of leaving them in a separate default `results/` directory. On macOS with the local profile, only `CLOOPS_CALL` is restricted to one CPU to avoid cLoops2 0.0.3 spawned-worker global-state errors; other tasks retain `hichip_loop_cpus`. The fix was checked using synthetic contacts and harmless configuration/resource tests, without restarting real-data analysis.

## Resolution parameters

- `base_resolution`: finest resolution used to create the initial `.cool` matrix.
- `resolutions`: resolutions stored in `.mcool`; every value must be a multiple of `base_resolution`.
- `analysis_resolution`: resolution used for coverage and insulation.
- `compartment_resolution`: resolution used for A/B compartment analysis.
- `loop_resolution`: resolution used for expected-cis, dots and FitHiChIP.

The finest configurable resolution is not necessarily biologically supportable. Select the final resolution using valid cis contacts, matrix sparsity, replicate concordance and feature recovery.

## Compute resources

`memory_per_cpu_gb` controls memory allocation per requested CPU and defaults to `4`. Each task receives `task.cpus * memory_per_cpu_gb` GB, while CPU count and time limit remain assigned by workload class:

| Class | CPU | Default memory | Time | Main processes |
|---|---:|---:|---:|---|
| `small` | 2 | 8 GB | 2 h | samplesheet validation, FastQC, contact QC and MultiQC |
| `medium` | 8 | 32 GB | 12 h | pairtools, BigWig, expected-cis, compartments, insulation/TADs and peak calling |
| `large` | 16 | 64 GB | 48 h | fastp, BWA/SAMtools, cooler and cooltools dots |

For a cluster with a different CPU-to-memory allocation ratio, change only `memory_per_cpu_gb` in the parameter file. The test profile intentionally overrides these production resources with 1 CPU and 1 GB per class.

## Outputs

```text
results/
├── 01_qc/                 FastQC, fastp, alignment and pairtools QC
├── 02_preprocessing/      Optional trimmed FASTQ and alignment BAM files
├── 03_pairs/              Selected UU/UR/RU pairs, optionally without deduplication
├── 04_matrices/           Base-resolution .cool and multi-resolution .mcool
├── 05_features/           Contact QC, expected, compartments, TADs and Hi-C loops
├── 06_hichip/             RPGC BigWig, MACS2 outputs, annotation; optional loop entry outputs
└── pipeline_info/         Validated samplesheet, MultiQC, trace, timeline, report and DAG
```

Alignment QC contains separate `SAMPLE.contacts.flagstat` and `SAMPLE.peaks.flagstat` files. The first summarizes the unfiltered BWA contact alignments, while the second summarizes the MAPQ-filtered, primary, non-duplicate alignments used for HiChIP peak calling. Both are included in MultiQC.

### Valid-pair retention

The contact branch selects `UU`, `UR` and `RU`, including both cis and trans contacts. `UR/RU` are successfully rescued, uniquely positioned contacts, not multi-mapping pairs. Ambiguous alignments and unmapped ends remain excluded.

`skip_pairtools_dedup` controls duplicate handling for the contact branch, independently of assay:

| Setting | Pairtools stages | Repeated contacts |
|---|---|---|
| `false` (default) | parse -> sort -> dedup -> select | Detected duplicates are excluded |
| `true` | parse -> sort -> select | Retained with their original UU/UR/RU types; no DD marking |

Append `--skip_pairtools_dedup true` to the main workflow command, or set `skip_pairtools_dedup: true` in the parameter YAML. This bypasses the entire `pairtools dedup` step, rather than using `--no-mark-dups`, which would only change labeling and would still remove duplicates. MAPQ, walks and selection settings are unchanged. UMI and linker processing are not added.

The task publishes `01_qc/pairtools/SAMPLE.dedup_policy.tsv` with `pairtools_dedup` set to `enabled` or `skipped`. `SAMPLE.dedup.stats` is generated only when deduplication runs; skipped duplicate counts are unknown, not a measured zero. Parse/final statistics and MultiQC remain available in both modes. Final statistics classify repeated unmarked records as UU/UR/RU and must not be interpreted as independent molecule counts.

After selection, `check_pairs_retention.py` compares every eligible type against `SAMPLE.dedup.stats` in default mode, or `SAMPLE.parse.stats` when deduplication is skipped. Any missing eligible records or unsupported retained type causes failure. `SAMPLE.retention.tsv` records expected and retained counts. A sample with no UR/RU is allowed.

Cooler aggregates all retained read pairs, including repeated coordinates when deduplication is skipped. Mustache and the COOL-input FitHiChIP routes use these rebuilt matrices. Native cLoops2 `pre -format pairs` still collapses identical coordinates internally; the switch does not disable cLoops2 deduplication. This limitation is printed in the loop-input QC log and recorded in the JSON report. The separate peak-calling BAM/BigWig branch is unchanged.

Previously generated UU-only or deduplicated pairs cannot recover discarded contacts after a code/configuration change. Deploy the updated checkout to the execution machine and rebuild pairs from upstream contact BAM/FASTQ, then matrices and loops, using a new output directory for each mode. `start_from: features` and the published-file loop entry point reuse existing matrices and cannot restore discarded contacts. Do not mix old matrices or stale QC files with a new deduplication policy. Results without deduplication are read-count-based exploratory results, not proven independent-molecule measurements.

For HiChIP samples, `save_peak_bigwig: true` generates `SAMPLE.peaks.RPGC.bw` from the filtered peak-calling BAM using deepTools `bamCoverage`. This follows the ChIPalign_NF RPGC/10-bp-bin convention, but intentionally omits `--extendReads`: HiChIP mates represent proximity contacts and must not be extended as contiguous ChIP-seq fragments. `bigwig_ignore_for_normalization` is optional; the hg38 example uses `chrX` to match the referenced ChIP-seq configuration.

MACS2 peak calling retains the ChIPalign_NF-style `SAMPLE_peaks.xls`, `SAMPLE_peaks.narrowPeak` and `SAMPLE_summits.bed` files under `06_hichip/callpeaks/SAMPLE.callpeaks/`, while `SAMPLE.peaks.bed` remains the canonical FitHiChIP input. For a supplied `peak_file`, the callpeaks directory instead records a normalized copy and its source path. When `run_peak_annotation` is enabled, HOMER `annotatePeaks.pl` writes `SAMPLE.peaks.anno.txt` and `SAMPLE.peaks.anno.stats` under `06_hichip/peak_annotation/`.

## Validation

Check pipeline parsing and DAG construction without executing analysis tools:

```bash
nextflow run main.nf -profile test -preview
```

Run the complete two-assay stub test:

```bash
nextflow run main.nf -profile test -stub-run
```

Validate a samplesheet independently:

```bash
python3 bin/validate_samplesheet.py \
  --input samplesheet.csv \
  --output /tmp/samplesheet.validated.csv \
  --base-dir "$PWD"
```

The stub test validates workflow scheduling and file routing only. It is not a biological validation or performance benchmark. Before production use, compare valid pairs, duplication, cis/trans ratios and contact-map sparsity against a previously validated dataset.

## Planned modules

Technical sequencing-run merging is implemented. The detailed roadmap is available in `docs/EXTENSIONS.md`. Priority additions are biological-replicate concordance and controlled replicate aggregation, depth-matched comparisons, differential contacts and loops, APA/pileup, regulatory annotation, and browser-ready `.hic` or HiGlass outputs.

## License

HiCLink_NF is released under the [GNU General Public License v3.0](LICENSE). You may use, modify and redistribute the workflow under the terms of the GPL-3.0 license. Distributed modified versions must retain the same license and make their corresponding source code available.
