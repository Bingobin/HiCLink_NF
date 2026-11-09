# HiChIP loop-calling and differential-loop subworkflow

This optional DSL2 subworkflow implements cLoops2, Mustache and Docker-based FitHiChIP analysis from existing pairs and cooler matrices. The FASTQ entry point (`main.nf`) publishes upstream inputs but does not call this subworkflow automatically. Its old direct FitHiChIP process and legacy parameters/profiles have been removed, so FitHiChIP runs only here. A separate entry point (`hichip_loops.nf`) supports reuse of published files without repeating upstream analysis. The FitHiChIP branch additionally needs peaks; see [FitHiChIP installation, union peaks and differential analysis](FITHICHIP_LOOPS.md).

`--run_loops false` disables the standalone entry point before input discovery or caller execution. `run_loops` defaults to true. `hichipLoopSettings(params)` includes this flag; directly supplied settings may use the Boolean `run_loops: false` to disable caller/preparation/differential branches while retaining input QC. No core workflow inputs are changed.

Repeated samplesheet IDs representing sequencing runs of one library are collapsed to one published library input. Assay, condition, replicate, control and peak file must agree across rows. FASTQ paths and generated `run_id` values are not used for loop analysis. Published pairs/matrices must already contain the merged runs; discovery cannot combine old per-run results.

Status: metadata/format tests, static linting, and previews were checked at initial setup. A subsequent user-run on macOS completed preparation but failed in cLoops2's spawned-worker logging path. The compatibility fix was verified with synthetic contacts and harmless resource/configuration jobs; real-data analysis was not restarted during debugging. No final biological result or performance benchmark is claimed.

## Standalone entry point

```bash
cp conf/hichip_loops_job.example.config conf/hichip_loops_job.local.config
# Edit the paths and optional hichip_loops_env prefix in the copied file.
nextflow run hichip_loops.nf \
  -c conf/hichip_loops_job.local.config -profile hichip_loops_local
```

Use one ordered config graph: `includeConfig 'hichip_loops.config'`, then job parameters or an included job file, then `includeConfig 'hichip_loops_runtime.config'`. The last include sees the final output setting and rebinds report, trace, timeline, and DAG paths. Separate `-c` files are parsed in separate parameter contexts; merely adding a report config last can still interpolate a missing/stale parameter. The machine-local launcher uses `hichip_loops.launch.config` to enforce the ordered graph while leaving the user's parameter file editable.

If `previous_results` is omitted, the directory above the samplesheet's parent is used. The expected layout is `RESULTS/pipeline_info/samplesheet.validated.csv`, `RESULTS/03_pairs/`, and `RESULTS/04_matrices/`. Explicit `previous_results` is recommended.

Only `ID` and `assay` are mandatory CSV columns. `condition`, `replicate`, and `control` are preserved; `control` contains the matched control **sample ID**, not a condition label. R1/R2 may point to unavailable server files: they are not accessed. FASTQ, BAM and reference FASTA are not required. cLoops2/Mustache need no peaks; FitHiChIP needs exact per-sample published peaks or an explicit absolute `peak_file`.

Both `hichip` and `chiapet` rows are selected as already processed contacts; `hic` rows are excluded. ChIA-PET input must already have pairs and a compatible matrix, plus peaks for FitHiChIP. This does not add raw-read/linker/UMI processing or establish protocol-specific statistical validation. Cross-assay differential comparisons are rejected. Each selected sample must have both pairs and a matrix, even if only one caller is enabled; this is the common subworkflow input contract. Exact ID matching, never substring matching, uses this priority:

1. Pairs: `03_pairs/ID.valid.pairs.gz`, then `ID.valid.pairs`, `ID.pairs.gz`, or `ID.pairs`.
2. Matrices: `04_matrices/ID.mcool`, then `ID.RESOLUTION.cool`, then `ID.cool`, then the closest finer `ID.BINSIZE.cool` whose bin size divides the requested resolution.
3. Actual HDF5 metadata is validated in a task. No upsampling or fractional coarsening is allowed. An incompatible preferred `.mcool` fails rather than silently selecting another product.

Blank `control` means individual calls only. Multiple treatments can reuse one control. Missing/self controls fail. For `CellLine.IP.Group` IDs, cross-cell-line/cross-IP contrasts fail; other naming conventions require the user to establish appropriate matching.

The current HMGA2 sheet selects five samples and these contrasts:

| Treatment | Control |
| --- | --- |
| A549.RAD21.SH2 | A549.RAD21.SCR |
| H1299.RAD21.SH23 | H1299.RAD21.SCR |

`A549.HMGA2.OE_HA` is called individually, not compared with RAD21 IP libraries.

## Interface for future main-workflow integration

The core reads no CSV and discovers no local directories:

| Input | Contract |
| --- | --- |
| `samples` | Queue channel: `tuple(meta, pairs_path, matrix_path)`; unique `meta.id`, `meta.assay` is `hichip` or `chiapet`; retain condition/replicate/control as needed |
| `comparisons` | Queue channel of maps: `[id: comparison_id, treatment: sample_id, control: sample_id]` |
| `cfg` | Ordinary immutable settings map, not a queue channel; use `hichipLoopSettings(params)` after loading the config, or supply the same settings directly |

Modules use tools from `PATH`; the local Conda prefix belongs to the optional profile, not the core. The shared helper is `bin/hichip_loops.py` in the parent project.

Example hook for the FASTQ branch of `main.nf` (documentation only, **not applied**):

```groovy
include { HICHIP_LOOPS; hichipLoopSettings } from './subworkflows/local/hichip_loops'

// In the parent config, include the base settings and hichip_loops_macos.config.
// Do not include the standalone runtime report rebinding in a full parent run.
loopSamples = PAIRTOOLS.out.pairs
    .filter { meta, pairs -> meta.assay == 'hichip' }
    .map { meta, pairs -> tuple(meta.id, meta, pairs) }
    .join(COOLER_MATRIX.out.matrices.map { meta, baseCool, multiCool -> tuple(meta.id, multiCool) })
    .map { id, meta, pairs, matrix -> tuple(meta, pairs, matrix) }

loopComparisons = sampleRecords
    .unique { meta, read1, read2 -> meta.id }
    .filter { meta, read1, read2 -> meta.assay == 'hichip' && meta.control }
    .map { meta, read1, read2 ->
        [id: "${meta.id}_vs_${meta.control}", treatment: meta.id, control: meta.control]
    }

HICHIP_LOOPS(loopSamples, loopComparisons, hichipLoopSettings(params))
```

The parent must supply every comparison sample in `samples`, unique comparison IDs, and appropriate biological matching. Collection holds metadata/path tuples, not contact arrays. In feature-restart mode, construct equivalent channels from published data, or use the separate discovery entry point.

Emitted channels: `input_qc`, `cloops_pre`, `cloops_calls`, `cloops_diff`, `mustache_cools`, `mustache_calls`, and `mustache_diff`. Disabled branches emit empty channels. Pre/call/difference outputs emit `tuple(meta_or_comparison, result_directory)`; Mustache cools emit `tuple(meta, working_cool, input_qc_report)`.

## Configuration

Load the base `conf/hichip_loops.config` through the ordered job config described above:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `hichip_loop_callers` | `cloops2,mustache` | Select from cloops2,mustache,fithichip; FitHiChIP setup is documented separately |
| `hichip_loop_diff` | `true` | Set false for individual calls only |
| `loop_resolution` | `10000` | Mustache resolution, inherited from the parent config |
| `min_loop_distance`, `max_loop_distance` | `20000`, `2000000` | PET distances for cLoops2; final bin-start distances for Mustache |
| `hichip_loop_chroms` | chr1-22,chrX | Cis only; exclude chrY/chrM |
| `hichip_loop_cpus` | `2` | CPUs per analysis task |
| `cloops2_eps`, `cloops2_minpts` | `5000,10000`, `10` | Initial clustering settings, not fixed-bin resolutions |
| `cloops2_diff_pcut`, `cloops2_diff_ma_fdr` | `0.01`, `0.1` | Bonferroni-adjusted Poisson threshold and MA background cutoff, not BH q |
| `mustache_pt`, `mustache_diff_pt`, `mustache_sparsity` | `0.05`, `0.05`, `0.88` | Native detection, differential, and sparsity thresholds |
| `loop_balance_max_iters` | `1000` | Working-copy balancing iteration limit |
| `hichip_loop_outdir` | `OUTDIR/06_hichip/loops` | Published loop outputs |

Mustache 1.3.3 silently adjusts some distances. This workflow instead requires `200*resolution <= max_loop_distance <= min(2000*resolution, 2 Mb)`. Start at 10 kb; run 5 kb sensitivity analysis separately if appropriate.

The `hichip_loops_local` profile leaves concurrency to Nextflow's resource-aware local scheduler: it does not set `queueSize` or `maxForks`. Nextflow detects the host CPU/memory capacity and admits ready tasks while their combined resource requests fit; concurrency is not fixed to a sample count. Task CPU/memory requests remain necessary estimates, not automatically measured requirements. This accounting does not track memory used by unrelated applications and is not a hard RAM limit. Analysis memory remains `task.cpus * memory_per_cpu_gb` GB, while FitHiChIP keeps its separate container resource limits; the Docker Desktop VM must also accommodate concurrently running containers. The profile limits numerical-library threads to one, keeps matplotlib caches task-local, and uses an existing environment via `hichip_loops_env`, or the current `PATH` if unset. It does not install an environment. Other executors can use the normal parent profiles with tools provisioned separately. Restart with `-resume` to apply configuration changes; an already-running session does not reload them.

The runtime config includes `hichip_loops_macos.config`, which also binds analysis CPU requests lazily: loading the base config before a local override must not freeze its default of two CPUs. Under the local profile on macOS, `CLOOPS_CALL` requests one CPU and consequently calls `cLoops2 callLoops -p 1`. In cLoops2 0.0.3, `logger` and `DBSCAN` are module globals initialized in the parent only; spawn workers re-import them as `None`. Both globals matter, so initializing only the logger would not repair the full defect. This policy avoids the faulty multiprocess path without forcing unsafe fork behavior or modifying site-packages. `CLOOPS_PRE`, `CLOOPS_DIFF`, cooler, and Mustache keep the final requested CPU count; Linux local runs also keep it. For future main-workflow integration, include the reusable resource/macOS policy separately, while leaving the parent's report locations under its own control.

Interfaces were checked against cLoops2 0.0.3, mustache-hic 1.3.3, cooler 0.12.1, Python 3.10, and Nextflow 26.04.1. The core Conda YAML is unchanged and does not provide the loop callers; supply a suitable environment.

## Safety, caching, and outputs

Original pairs/matrices remain unchanged. Mustache copies the selected layer or coarsens into a new cool. Converged weights are reused; missing/unconverged weights are balanced only in that copy using `ignore-diags=2`, `mad-max=5`, `min-nnz=10`, `min-count=0`, `tol=1e-5`. Failed convergence stops the task, without a raw-count fallback. Read-only inspection after the user's partial run confirmed that the original H1299-SCR 10 kb layer remains unconverged, while the prepared working cool has converged weights. Existing prepared data and caches are retained.

cLoops2 pre metadata embeds absolute paths. Call/difference tasks stage **copies** and run `cLoops2 update` on those copies, never upstream cached directories. Keep the work directory for `-resume`; published pre datasets moved elsewhere need `cLoops2 update` for standalone reuse. No loop-based PET filtering or automatic downsampling is applied.

Outputs are `input_qc/`, `cloops2/{pre,call,diff}/`, `mustache/{cool,call,diff}/`, and `pipeline_info/`. Discovery is deliberately uncached to recheck file stats. Its signature guard refuses the same output directory after changes to settings, samplesheet contents, helper code, or input sizes/mtimes. Use a **new output directory** for changed parameters/inputs; use `-resume` for unchanged retries. Published files are not overwritten. Signatures are metadata-based, not full matrix/pairs checksums. When bypassing discovery, the integrating parent is responsible for the equivalent output-directory policy.

The standalone runtime config now writes all four execution reports into `hichip_loop_outdir/pipeline_info/`. Old reports already written under `results/pipeline_info/` are historical outputs and are not deleted or moved by the fix. The discovery helper and scientific settings were not changed by this debug patch, and the existing HMGA2 input signature still passes. `-resume` can reuse successful matching tasks; changing CPU counts or commands may invalidate the affected tasks, so not every previous task is guaranteed to be cached.

For cLoops2, interpret `significant` together with `scaledFc` in `*_dloops.txt`: positive means treatment-enriched, negative means control-enriched. Native specific-output names refer to staging roles `treatment_pre` and `control_pre`; `comparison.tsv` records their sample IDs.

For differential Mustache, f1=treatment and f2=control: `.diffloop1` is treatment-enriched, `.diffloop2` control-enriched. Native and distance-filtered `*.treatment_enriched.candidates.tsv` / `*.control_enriched.candidates.tsv` are retained. The native FDR is not relabeled as a replicate-aware differential q-value. Header-only Mustache output is valid; missing/invalid files fail. cLoops2 stops difference testing if either library has no candidate loops. Mustache's import issue is handled by executing the package's `diff_mustache.py` directly.

## Interpretation and validation limits

Current HMGA2 contrasts have one library per group. These methods produce **exploratory candidates**, not biological-replicate-aware inference. Do not split reads into pseudo-biological replicates, equate group-exclusive detection with differential significance, compare different IPs/cell lines directly, or interpret every HiChIP change as structural rewiring. Depth, RAD21 occupancy, IP efficiency, and accessibility affect signal. Add biological repeats, depth-matching sensitivity analyses, APA/local maps, and ChIP/ATAC/RNA evidence before mechanistic conclusions. Mustache was developed for Hi-C/Micro-C; its model does not automatically remove HiChIP enrichment bias.

Checks without executing analysis tasks:

```bash
nextflow lint hichip_loops.nf subworkflows/local/hichip_loops.nf \
  modules/local/hichip_loops.nf conf/hichip_loops*.config tests/loop_runtime.nf
python -m unittest discover -s tests -p 'test_*.py' -v
nextflow run tests/loop_channels.nf  # channel routing only; no processes or data analysis
nextflow run hichip_loops.nf -c conf/hichip_loops_job.local.config -preview
```

`tests/loop_runtime.nf` uses two harmless shell tasks named like the caller processes to verify CPU routing (macOS: cLoops2=1, Mustache=requested CPUs), the environment prefix, and the requested memory; it never calls either analysis tool. Run it with the ordered job config and an isolated temporary work/output directory, not the production output directory.

Optional native-tool smoke test (synthetic contacts only):

```bash
python tests/smoke_cloops2.py --executable /path/to/env/bin/cLoops2 --check-parallel-failure
```

On the observed macOS environment this reproduced the `-p 2` logger error and generated a valid loop table with `-p 1` from 840 synthetic contacts. No stubs are supplied for the real entry point; do not use `-stub-run` instead of a no-task preview. Complete real-data downstream execution, RAM peaks, optimal parameters, and biological recovery remain unvalidated after this fix.

References: [cLoops2](https://github.com/YaqiangCao/cLoops2), [Mustache](https://github.com/ay-lab/mustache), [cooler](https://cooler.readthedocs.io/en/latest/cli.html), and [Nextflow processes](https://docs.seqera.io/nextflow/process/).
