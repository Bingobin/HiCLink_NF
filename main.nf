#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_SAMPLESHEET; FASTQC_RAW; FASTP; BWA_ALIGN; PEAKS_BIGWIG; PAIRTOOLS; COOLER_MATRIX } from './modules/local/preprocessing'
include { CONTACT_QC; EXPECTED_CIS; COMPARTMENTS; INSULATION_TADS; COOLTOOLS_DOTS; HICHIP_PEAKS; HICHIP_PEAK_ANNOTATION; FITHICHIP } from './modules/local/features'
include { MULTIQC } from './modules/local/reporting'


def intList(value, String name) {
    def values = value instanceof Collection ? value : value.toString().split(',')
    try {
        return values.collect { it.toString().trim() }
                     .findAll { it }
                     .collect { it as Integer }
                     .unique()
                     .sort()
    } catch (Exception ignored) {
        error "Parameter --${name} must be a comma-separated list of integers; received: ${value}"
    }
}


workflow {
    int memoryPerCpuGb
    try {
        memoryPerCpuGb = params.memory_per_cpu_gb as Integer
    } catch (Exception ignored) {
        error "Parameter --memory_per_cpu_gb must be a positive integer; received: ${params.memory_per_cpu_gb}"
    }
    if (memoryPerCpuGb <= 0) {
        error "Parameter --memory_per_cpu_gb must be a positive integer; received: ${params.memory_per_cpu_gb}"
    }

    int macs2Extsize
    try {
        macs2Extsize = params.macs2_extsize as Integer
    } catch (Exception ignored) {
        error "Parameter --macs2_extsize must be a positive integer; received: ${params.macs2_extsize}"
    }
    if (macs2Extsize <= 0) {
        error "Parameter --macs2_extsize must be a positive integer; received: ${params.macs2_extsize}"
    }

    log.info """
    HiCLink_NF ${workflow.manifest.version}
    Unified Hi-C and HiChIP analysis
    =================================
    Samplesheet : ${params.input}
    Assay mode  : determined per sample
    Output      : ${params.outdir}
    Genome      : ${params.genome}
    Memory/CPU  : ${memoryPerCpuGb} GB
    """.stripIndent()

    if (!params.input) {
        error 'A samplesheet is required. Set --input <samplesheet.csv>.'
    }
    if (!params.fasta) {
        error 'A BWA-indexed reference FASTA is required. Set --fasta <reference.fa>.'
    }
    if (!params.chrom_sizes) {
        error 'A chromosome sizes file is required. Set --chrom_sizes <genome.chrom.sizes>.'
    }

    def resolutions = intList(params.resolutions, 'resolutions')
    def tadWindows  = intList(params.tad_windows, 'tad_windows')
    def requiredResolutions = [params.base_resolution as Integer, params.analysis_resolution as Integer]
    if (params.run_compartments.toString().toBoolean()) {
        requiredResolutions << (params.compartment_resolution as Integer)
    }
    if (params.run_hic_loops.toString().toBoolean() || params.run_fithichip.toString().toBoolean()) {
        requiredResolutions << (params.loop_resolution as Integer)
    }
    def missingResolutions = requiredResolutions.unique().findAll { !resolutions.contains(it) }
    if (missingResolutions) {
        error "The following requested resolutions are absent from --resolutions: ${missingResolutions.join(',')}"
    }
    def invalidMultiples = resolutions.findAll { it % (params.base_resolution as Integer) != 0 }
    if (invalidMultiples) {
        error "Every zoom resolution must be a multiple of --base_resolution; invalid: ${invalidMultiples.join(',')}"
    }
    if (!tadWindows) {
        error '--tad_windows must contain at least one window size.'
    }
    if (params.run_compartments.toString().toBoolean() && !params.phasing_track) {
        log.warn 'Compartment analysis was requested without --phasing_track; COMPARTMENTS will be skipped.'
    }

    inputSheet = file(params.input, checkIfExists: true)
    referenceFasta = file(params.fasta, checkIfExists: true)
    chromSizes = file(params.chrom_sizes, checkIfExists: true)
    emptyPeak = file("${projectDir}/assets/empty.bed", checkIfExists: true)
    multiqcConfig = file(params.multiqc_config, checkIfExists: true)

    // Stage the classic BWA index beside the FASTA so that local and container
    // execution behave identically.
    bwaIndexes = ['.amb', '.ann', '.bwt', '.pac', '.sa'].collect { suffix ->
        file("${params.fasta}${suffix}", checkIfExists: true)
    }

    VALIDATE_SAMPLESHEET(
        Channel.value(inputSheet),
        Channel.value(inputSheet.parent.toString())
    )

    samples = VALIDATE_SAMPLESHEET.out.csv
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def peakFile = row.peak_file?.trim() ?: ''
            def meta = [
                id            : row.ID,
                assay         : row.assay,
                condition     : row.condition,
                replicate     : row.replicate,
                control       : row.control,
                peak_file     : peakFile,
                peak_provided : peakFile as Boolean
            ]
            tuple(
                meta,
                file(row.R1, checkIfExists: true),
                file(row.R2, checkIfExists: true)
            )
        }

    FASTQC_RAW(samples)
    FASTP(samples)

    alignmentInput = FASTP.out.reads.map { meta, read1, read2 ->
        tuple(meta, read1, read2, referenceFasta, bwaIndexes)
    }
    BWA_ALIGN(alignmentInput)

    hichipPeakBams = BWA_ALIGN.out.peak_bam
        .filter { meta, bam, bai -> meta.assay == 'hichip' }
    PEAKS_BIGWIG(hichipPeakBams)

    pairtoolsInput = BWA_ALIGN.out.contacts.map { meta, bam ->
        tuple(meta, bam, chromSizes)
    }
    PAIRTOOLS(pairtoolsInput)

    coolerInput = PAIRTOOLS.out.pairs.map { meta, pairs ->
        tuple(meta, pairs, chromSizes)
    }
    COOLER_MATRIX(coolerInput)

    matrices = COOLER_MATRIX.out.matrices
    matrixUris = matrices.map { meta, baseCool, multiCool -> tuple(meta, multiCool) }
    globalFeatureMatrices = matrixUris.filter { meta, multiCool ->
        meta.assay == 'hic' || params.global_features_all_assays.toString().toBoolean()
    }
    expectedInput = matrixUris.filter { meta, multiCool ->
        meta.assay == 'hic' ||
            params.global_features_all_assays.toString().toBoolean() ||
            params.dots_all_assays.toString().toBoolean()
    }

    CONTACT_QC(matrixUris)
    EXPECTED_CIS(expectedInput)

    phasingTrack = params.phasing_track
        ? file(params.phasing_track, checkIfExists: true)
        : emptyPeak
    compartmentInput = globalFeatureMatrices.map { meta, multiCool ->
        tuple(meta, multiCool, phasingTrack)
    }
    COMPARTMENTS(compartmentInput)
    INSULATION_TADS(globalFeatureMatrices)

    hicLoopInput = EXPECTED_CIS.out.expected
        .filter { meta, multiCool, expected ->
            meta.assay == 'hic' || params.dots_all_assays.toString().toBoolean()
        }
    COOLTOOLS_DOTS(hicLoopInput)

    hichipPeakInput = hichipPeakBams.map { meta, bam, bai ->
            def peak = meta.peak_provided
                ? file(meta.peak_file, checkIfExists: true)
                : emptyPeak
            tuple(meta, bam, bai, peak)
    }
    HICHIP_PEAKS(hichipPeakInput)
    HICHIP_PEAK_ANNOTATION(HICHIP_PEAKS.out.peaks)

    keyedPeaks = HICHIP_PEAKS.out.peaks.map { meta, peaks -> tuple(meta.id, meta, peaks) }
    keyedMatrices = matrices
        .filter { meta, baseCool, multiCool -> meta.assay == 'hichip' }
        .map { meta, baseCool, multiCool -> tuple(meta.id, multiCool) }
    fithichipInput = keyedPeaks.join(keyedMatrices).map { id, meta, peaks, multiCool ->
        tuple(meta, peaks, multiCool, chromSizes)
    }
    FITHICHIP(fithichipInput)

    multiqcInputs = FASTQC_RAW.out.zips
        .map { meta, files -> files }
        .flatten()
        .mix(FASTP.out.json.map { meta, json -> json })
        .mix(PAIRTOOLS.out.stats.map { meta, parseStats, dedupStats, pairStats -> [parseStats, dedupStats, pairStats] }.flatten())
        .mix(BWA_ALIGN.out.flagstats.map { meta, contactsFlagstat, peaksFlagstat ->
            [contactsFlagstat, peaksFlagstat]
        }.flatten())

    MULTIQC(multiqcInputs.collect(), Channel.value(multiqcConfig))

    workflow.onComplete = {
        def headline = workflow.success
            ? "Done! See results --> ${params.outdir}"
            : "Pipeline failed! Check .nextflow.log for details."

        def summary = """
        ${headline}

        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Status      : ${workflow.success ? 'SUCCEEDED' : 'FAILED'}
        Exit status : ${workflow.exitStatus}
        """.stripIndent().trim()

        if (!workflow.success && workflow.errorMessage) {
            summary += "\nError       : ${workflow.errorMessage}"
        }

        log.info "\n${summary}\n"
    }
}
