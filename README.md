# HiCLink_NF

`HiCLink_NF` is a Nextflow DSL2 pipeline for unified Hi-C and HiChIP analysis.

The original `hic_cool.nf` and `hic_proc.nf` scripts are retained as legacy entry points. The new workflow starts from `main.nf`.

## Workflow

```text
FASTQ
  ├─ FastQC and fastp
  ├─ BWA-MEM alignment (-SP)
  ├─ pairtools parse, sort, deduplication, filtering and QC
  ├─ balanced single- and multi-resolution cooler matrices
  ├─ shared: contact coverage and matrix sparsity QC
  ├─ Hi-C: expected-cis, compartments, insulation/TADs and cooltools dots
  └─ HiChIP: RPGC BigWig, supplied or MACS2 peaks, HOMER annotation and FitHiChIP loops
```

A single samplesheet can contain both `hic` and `hichip` libraries. The pipeline retains `condition`, `replicate` and `control` metadata for future replicate-aware and differential modules.

HiChIP matrices contain target-enrichment bias. Compartment and TAD analysis therefore runs on Hi-C samples by default. Use `--global_features_all_assays true` only when there is a specific methodological reason to calculate these features for HiChIP, and interpret them separately from unbiased Hi-C maps.

## Samplesheet

Use `assets/samplesheet.example.csv` as the template.

| Column | Required | Description |
|---|---:|---|
| `ID` | Yes | Unique sample ID containing letters, numbers, `.`, `_` or `-` |
| `R1`, `R2` | Yes | Paired-end FASTQ files; relative paths are resolved from the samplesheet directory |
| `assay` | Yes | `hic` or `hichip` |
| `condition` | No | Biological condition; defaults to the sample ID |
| `replicate` | No | Biological replicate number; defaults to 1 |
| `control` | No | Sample ID of a matched control, reserved for differential analysis |
| `peak_file` | HiChIP only | External ChIP-seq or consensus peaks; MACS2 is used when this field is empty |

Example:

```csv
ID,R1,R2,assay,condition,replicate,control,peak_file
HIC_CTRL_R1,/data/hic_ctrl_R1.fastq.gz,/data/hic_ctrl_R2.fastq.gz,hic,CTRL,1,,
HICHIP_TREAT_R1,/data/hichip_R1.fastq.gz,/data/hichip_R2.fastq.gz,hichip,TREAT,1,HICHIP_CTRL_R1,/data/H3K27ac.peaks.bed
```

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

FitHiChIP is kept outside the core Conda environment. Either provide a local `FitHiChIP_HiCPro.sh` through `--fithichip_command`, or use the supplied Docker/Apptainer FitHiChIP profile.

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

Run on SLURM and use the FitHiChIP Apptainer image:

```bash
nextflow run main.nf \
  -profile conda,slurm,apptainer_fithichip \
  -params-file conf/params.local.yml \
  --input samplesheet.csv
```

For a Hi-C-only project, set `--run_fithichip false`. For a HiChIP sample without `peak_file`, set `--macs2_genome_size` to `hs`, `mm` or an effective genome size. The workflow converts primary HiChIP alignments to BED before MACS2 peak calling because proximity-ligation pairs are not ordinary short-insert `BAMPE` fragments. MACS2 runs with `--nomodel`; `macs2_extsize` controls the fixed extension length and defaults to 147 bp.

### Restart feature analysis from published results

If `work/` was deleted, start a new feature run using the published results rather than relying on `-resume`:

```bash
nextflow run main.nf \
  -profile standard,apptainer_fithichip \
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

Missing matrices or required peak inputs produce an error naming the sample. For matrix features only, set `--run_hichip_peaks false --run_peak_annotation false --run_fithichip false`. When running FitHiChIP or peak annotation, enable `run_hichip_peaks` even if peaks are already available: this step imports the peak file instead of calling MACS2. Structural HiChIP features still require `--global_features_all_assays true`; the default Hi-C/HiChIP routing is unchanged.

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
| `large` | 16 | 64 GB | 48 h | fastp, BWA/SAMtools, cooler, cooltools dots and FitHiChIP |

For a cluster with a different CPU-to-memory allocation ratio, change only `memory_per_cpu_gb` in the parameter file. The test profile intentionally overrides these production resources with 1 CPU and 1 GB per class.

## Outputs

```text
results/
├── 01_qc/                 FastQC, fastp, alignment and pairtools QC
├── 02_preprocessing/      Optional trimmed FASTQ and alignment BAM files
├── 03_pairs/              Deduplicated UU/UR/RU valid pairs
├── 04_matrices/           Base-resolution .cool and multi-resolution .mcool
├── 05_features/           Contact QC, expected, compartments, TADs and Hi-C loops
├── 06_hichip/             RPGC BigWig, MACS2 outputs, peak annotation and FitHiChIP loops
└── pipeline_info/         Validated samplesheet, MultiQC, trace, timeline, report and DAG
```

Alignment QC contains separate `SAMPLE.contacts.flagstat` and `SAMPLE.peaks.flagstat` files. The first summarizes the unfiltered BWA contact alignments, while the second summarizes the MAPQ-filtered, primary, non-duplicate alignments used for HiChIP peak calling. Both are included in MultiQC.

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

The detailed roadmap is available in `docs/EXTENSIONS.md`. Priority additions are lane and replicate merging, replicate concordance, depth-matched comparisons, differential contacts and loops, APA/pileup, regulatory annotation, and browser-ready `.hic` or HiGlass outputs.

## License

HiCLink_NF is released under the [GNU General Public License v3.0](LICENSE). You may use, modify and redistribute the workflow under the terms of the GPL-3.0 license. Distributed modified versions must retain the same license and make their corresponding source code available.
