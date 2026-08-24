process CONTACT_QC {
    tag "${meta.id}"
    label 'small'

    publishDir "${params.outdir}/05_features/contact_qc", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(multi_cool)

    output:
    tuple val(meta), path("${meta.id}.contact_qc.tsv"), path("${meta.id}.coverage.tsv"), emit: qc

    when:
    params.run_contact_qc.toString().toBoolean()

    script:
    def uri = "${multi_cool}::/resolutions/${params.analysis_resolution}"
    """
    cool_sparsity.py '${uri}' --output ${meta.id}.contact_qc.tsv
    cooltools coverage \\
        --nproc ${task.cpus} \\
        -o ${meta.id}.coverage.tsv \\
        '${uri}'
    """

    stub:
    """
    touch ${meta.id}.contact_qc.tsv ${meta.id}.coverage.tsv
    """
}


process EXPECTED_CIS {
    tag "${meta.id}"
    label 'medium'

    publishDir "${params.outdir}/05_features/expected", pattern: '*.expected_cis.tsv', mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(multi_cool)

    output:
    tuple val(meta), path(multi_cool), path("${meta.id}.expected_cis.tsv"), emit: expected

    when:
    params.run_expected.toString().toBoolean() || params.run_hic_loops.toString().toBoolean()

    script:
    def uri = "${multi_cool}::/resolutions/${params.loop_resolution}"
    """
    cooltools expected-cis \\
        --nproc ${task.cpus} \\
        -o ${meta.id}.expected_cis.tsv \\
        '${uri}'
    """

    stub:
    """
    touch ${meta.id}.expected_cis.tsv
    """
}


process COMPARTMENTS {
    tag "${meta.id}"
    label 'medium'

    publishDir "${params.outdir}/05_features/compartments", pattern: '*.compartments*', mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(multi_cool), path(phasing_track)

    output:
    tuple val(meta), path("${meta.id}.compartments*"), emit: compartments

    when:
    params.run_compartments.toString().toBoolean() && params.phasing_track

    script:
    def uri = "${multi_cool}::/resolutions/${params.compartment_resolution}"
    def track = params.phasing_track_column
        ? "${phasing_track}::${params.phasing_track_column}"
        : phasing_track.toString()
    """
    cooltools eigs-cis \\
        --n-eigs ${params.compartment_n_eigs} \\
        --phasing-track '${track}' \\
        --out-prefix ${meta.id}.compartments \\
        '${uri}'
    """

    stub:
    """
    touch ${meta.id}.compartments.cis.vecs.tsv
    """
}


process INSULATION_TADS {
    tag "${meta.id}"
    label 'medium'

    publishDir "${params.outdir}/05_features/tads", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(multi_cool)

    output:
    tuple val(meta), path("${meta.id}.insulation.tsv"), path("${meta.id}.TAD.*bp.bed"), emit: tads

    when:
    params.run_tads.toString().toBoolean()

    script:
    def uri = "${multi_cool}::/resolutions/${params.analysis_resolution}"
    def windows = params.tad_windows instanceof Collection
        ? params.tad_windows.join(' ')
        : params.tad_windows.toString().split(',').collect { it.trim() }.join(' ')
    """
    cooltools insulation \\
        --nproc ${task.cpus} \\
        --threshold ${params.insulation_threshold} \\
        -o ${meta.id}.insulation.tsv \\
        '${uri}' \\
        ${windows}

    for window in ${windows}; do
        insulation_to_TAD.py \\
            ${meta.id}.insulation.tsv \\
            --window \${window} \\
            --output ${meta.id}.TAD.\${window}bp.bed
    done
    """

    stub:
    """
    touch ${meta.id}.insulation.tsv ${meta.id}.TAD.stub_bp.bed
    """
}


process COOLTOOLS_DOTS {
    tag "${meta.id}"
    label 'large'

    publishDir "${params.outdir}/05_features/loops/cooltools", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(multi_cool), path(expected)

    output:
    tuple val(meta), path("${meta.id}.dots.tsv"), emit: loops

    when:
    params.run_hic_loops.toString().toBoolean()

    script:
    def uri = "${multi_cool}::/resolutions/${params.loop_resolution}"
    """
    cooltools dots \\
        --nproc ${task.cpus} \\
        --max-loci-separation ${params.max_loop_distance} \\
        -o ${meta.id}.dots.tsv \\
        '${uri}' \\
        ${expected}
    """

    stub:
    """
    touch ${meta.id}.dots.tsv
    """
}


process HICHIP_PEAKS {
    tag "${meta.id}"
    label 'medium'

    publishDir "${params.outdir}/06_hichip/peaks", pattern: '*.peaks.bed', mode: params.publish_mode, overwrite: true
    publishDir "${params.outdir}/06_hichip/callpeaks", pattern: '*.callpeaks', mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(bam), path(bai), path(input_peak)

    output:
    tuple val(meta), path("${meta.id}.peaks.bed"), emit: peaks
    tuple val(meta), path("${meta.id}.callpeaks"), emit: callpeak_results

    when:
    params.run_hichip_peaks.toString().toBoolean()

    script:
    if (meta.peak_provided) {
        return """
        mkdir -p ${meta.id}.callpeaks
        if [[ '${input_peak}' == *.gz ]]; then
            gzip -cd ${input_peak} | awk 'BEGIN{OFS="\\t"} !/^#/ && NF >= 3 {print}' > ${meta.id}.peaks.bed
        else
            awk 'BEGIN{OFS="\\t"} !/^#/ && NF >= 3 {print}' ${input_peak} > ${meta.id}.peaks.bed
        fi
        test -s ${meta.id}.peaks.bed
        cp ${meta.id}.peaks.bed ${meta.id}.callpeaks/${meta.id}.external.peaks.bed
        printf 'source\\t%s\\n' '${meta.peak_file}' > ${meta.id}.callpeaks/${meta.id}.peak_source.tsv
        """
    }
    if (!params.macs2_genome_size) {
        error "HiChIP sample ${meta.id} has no peak_file; set --macs2_genome_size (for example hs, mm, or an effective genome size)."
    }
    """
    set -o pipefail

    # Hi-C/HiChIP alignments are proximity-ligation contacts, not ordinary
    # short-insert paired-end fragments. Convert both primary read alignments
    # to BED so MACS2 does not discard them in BAMPE mode.
    bedtools bamtobed -i ${bam} \\
        | awk '{ seen = 1; print } END { if (!seen) { print "ERROR: no primary alignments were available for MACS2 peak calling" > "/dev/stderr"; exit 42 } }' \\
        | pigz -p ${task.cpus} > ${meta.id}.primary.bed.gz

    macs2 callpeak \\
        --treatment ${meta.id}.primary.bed.gz \\
        --format BED \\
        --gsize ${params.macs2_genome_size} \\
        --name ${meta.id} \\
        --qvalue ${params.macs2_qvalue} \\
        --nomodel \\
        --extsize ${params.macs2_extsize} \\
        --SPMR \\
        --keep-dup all \\
        ${params.macs2_extra_args}

    cp ${meta.id}_peaks.narrowPeak ${meta.id}.peaks.bed
    mkdir -p ${meta.id}.callpeaks
    mv \\
        ${meta.id}_peaks.xls \\
        ${meta.id}_peaks.narrowPeak \\
        ${meta.id}_summits.bed \\
        ${meta.id}.callpeaks/
    rm -f ${meta.id}.primary.bed.gz
    """

    stub:
    """
    mkdir -p ${meta.id}.callpeaks
    touch ${meta.id}.peaks.bed
    touch \\
        ${meta.id}.callpeaks/${meta.id}_peaks.xls \\
        ${meta.id}.callpeaks/${meta.id}_peaks.narrowPeak \\
        ${meta.id}.callpeaks/${meta.id}_summits.bed
    """
}


process HICHIP_PEAK_ANNOTATION {
    tag "${meta.id}"
    label 'small'

    publishDir "${params.outdir}/06_hichip/peak_annotation", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(peaks)

    output:
    tuple val(meta), path("${meta.id}.peaks.anno.txt"), path("${meta.id}.peaks.anno.stats"), emit: annotation

    when:
    params.run_peak_annotation.toString().toBoolean() && meta.assay == 'hichip'

    script:
    """
    annotatePeaks.pl \\
        ${peaks} \\
        ${params.genome} \\
        -annStats ${meta.id}.peaks.anno.stats \\
        > ${meta.id}.peaks.anno.txt
    """

    stub:
    """
    touch ${meta.id}.peaks.anno.txt ${meta.id}.peaks.anno.stats
    """
}


process FITHICHIP {
    tag "${meta.id}"
    label 'xlarge'

    publishDir "${params.outdir}/06_hichip/fithichip", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(peaks), path(multi_cool), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.id}.fithichip"), emit: loops
    tuple val(meta), path("${meta.id}.fithichip.config.txt"), emit: config

    when:
    params.run_fithichip.toString().toBoolean()

    script:
    def uri = "${multi_cool}::/resolutions/${params.loop_resolution}"
    """
    mkdir -p ${meta.id}.fithichip

    cat > ${meta.id}.fithichip.config.txt <<'FITHICHIP_CONFIG'
ValidPairs=
Interval=
Matrix=
Bed=
HIC=
COOL=${uri}
PeakFile=${peaks}
OutDir=${meta.id}.fithichip
ChrSizeFile=${chrom_sizes}
CircularGenome=0
IntType=${params.fithichip_interaction_type}
BINSIZE=${params.loop_resolution}
LowDistThr=${params.min_loop_distance}
UppDistThr=${params.max_loop_distance}
UseP2PBackgrnd=${params.fithichip_use_p2p_background}
BiasType=${params.fithichip_bias_type}
MergeInt=${params.fithichip_merge_interactions}
QVALUE=${params.fithichip_qvalue}
PREFIX=${meta.id}
OverWrite=1
FITHICHIP_CONFIG

    ${params.fithichip_command} -C ${meta.id}.fithichip.config.txt
    """

    stub:
    """
    mkdir -p ${meta.id}.fithichip
    touch ${meta.id}.fithichip/stub.loop.txt ${meta.id}.fithichip.config.txt
    """
}
