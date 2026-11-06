#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HICHIP_LOOPS; hichipLoopSettings } from './subworkflows/local/hichip_loops'
include { DISCOVER_LOOP_INPUTS } from './modules/local/hichip_loops'


workflow {
    if (!(params.run_loops.toString().toLowerCase() in ['true', 'false'])) {
        error 'Parameter --run_loops must be true or false.'
    }
    if (params.help) {
        log.info '''
HiCLink_NF: standalone entry point for the reusable HiChIP loop subworkflow

nextflow run hichip_loops.nf -c conf/hichip_loops.config \
    --input /analysis/pipeline_info/samplesheet.validated.csv \
    --previous_results /analysis --outdir results_loops

Use -profile hichip_loops_local --hichip_loops_env /path/to/conda/envs/hiclink
to use an existing local environment. See docs/HICHIP_LOOPS.md for integration.
No FASTQ, BAM or FASTA is required. FitHiChIP additionally requires peaks for each sample.
'''
    } else if (!params.run_loops.toString().toBoolean()) {
        log.info 'Loop calling disabled (--run_loops false); no loop tasks or input discovery will run.'
    } else {
        if (!params.input) {
            error 'A samplesheet is required: --input <samplesheet.validated.csv>'
        }
        sheet = file(params.input, checkIfExists: true)
        publishedRoot = params.previous_results
            ? file(params.previous_results, checkIfExists: true)
            : sheet.parent.parent
        if (!publishedRoot.isDirectory()) {
            error "Published results root is not a directory: ${publishedRoot}"
        }
        cfg = hichipLoopSettings(params)
        helper = file("${projectDir}/bin/hichip_loops.py", checkIfExists: true)
        log.info "HiChIP loop input: ${sheet}; published files: ${publishedRoot}; output: ${cfg.outdir}"
        log.warn 'Single-library differential loops are exploratory; no biological replicate variance is estimated.'

        DISCOVER_LOOP_INPUTS(sheet, publishedRoot.toString(), cfg, helper)
        samples = DISCOVER_LOOP_INPUTS.out.samples.splitCsv(header: true, sep: '\t')
            .map { row ->
                def meta = [id: row.ID, assay: row.assay, condition: row.condition,
                            replicate: row.replicate, control: row.control, peak_file: row.peak_file ?: '']
                tuple(meta, file(row.pairs, checkIfExists: true), file(row.matrix, checkIfExists: true))
            }
        comparisons = DISCOVER_LOOP_INPUTS.out.comparisons.splitCsv(header: true, sep: '\t')
            .map { row -> [id: row.id, treatment: row.treatment, control: row.control] }
        HICHIP_LOOPS(samples, comparisons, cfg)
    }
}
