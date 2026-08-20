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
  └─ HiChIP: supplied or MACS2 peaks and FitHiChIP loops
```

A single samplesheet can contain both `hic` and `hichip` libraries. The pipeline retains `condition`, `replicate` and `control` metadata for future replicate-aware and differential modules.

HiChIP matrices contain target-enrichment bias. Compartment and TAD analysis therefore runs on Hi-C samples by default. Use `--global_features_all_assays true` only when there is a specific methodological reason to calculate these features for HiChIP, and interpret them separately from unbiased Hi-C maps.

## Samplesheet

Use `assets/samplesheet.example.csv` as the template.

| Column | Required | Description |
|---|---:|---|
| `sample` | Yes | Unique sample ID containing letters, numbers, `.`, `_` or `-` |
| `assay` | Yes | `hic` or `hichip` |
| `fastq_1`, `fastq_2` | Yes | Paired-end FASTQ files; relative paths are resolved from the samplesheet directory |
| `condition` | No | Biological condition; defaults to the sample ID |
| `replicate` | No | Biological replicate number; defaults to 1 |
| `control` | No | Sample ID of a matched control, reserved for differential analysis |
| `peak_file` | HiChIP only | External ChIP-seq or consensus peaks; MACS2 is used when this field is empty |

Example:

```csv
sample,assay,fastq_1,fastq_2,condition,replicate,control,peak_file
HIC_CTRL_R1,hic,/data/hic_ctrl_R1.fastq.gz,/data/hic_ctrl_R2.fastq.gz,CTRL,1,,
HICHIP_TREAT_R1,hichip,/data/hichip_R1.fastq.gz,/data/hichip_R2.fastq.gz,TREAT,1,HICHIP_CTRL_R1,/data/H3K27ac.peaks.bed
```

## Reference files

Required reference files:

- A reference FASTA with classic BWA index files beside it: `.amb`, `.ann`, `.bwt`, `.pac` and `.sa`.
- A two-column chromosome sizes file without a header.
- An optional four-column, header-containing phasing track for compartment analysis. Its resolution must match `compartment_resolution`, and `phasing_track_column` identifies its signal column, for example `GC`.

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

For a Hi-C-only project, set `--run_fithichip false`. For a HiChIP sample without `peak_file`, set `--macs2_genome_size` to `hs`, `mm` or an effective genome size.

## Resolution parameters

- `base_resolution`: finest resolution used to create the initial `.cool` matrix.
- `resolutions`: resolutions stored in `.mcool`; every value must be a multiple of `base_resolution`.
- `analysis_resolution`: resolution used for coverage and insulation.
- `compartment_resolution`: resolution used for A/B compartment analysis.
- `loop_resolution`: resolution used for expected-cis, dots and FitHiChIP.

The finest configurable resolution is not necessarily biologically supportable. Select the final resolution using valid cis contacts, matrix sparsity, replicate concordance and feature recovery.

## Outputs

```text
results/
├── 01_qc/                 FastQC, fastp, alignment and pairtools QC
├── 02_preprocessing/      Optional trimmed FASTQ and alignment BAM files
├── 03_pairs/              Deduplicated UU valid pairs
├── 04_matrices/           Base-resolution .cool and multi-resolution .mcool
├── 05_features/           Contact QC, expected, compartments, TADs and Hi-C loops
├── 06_hichip/             Peaks and FitHiChIP loops
└── pipeline_info/         Validated samplesheet, MultiQC, trace, timeline, report and DAG
```

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
