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
    cooltools coverage \
        --nproc ${task.cpus} \
        -o ${meta.id}.coverage.tsv \
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
    cooltools expected-cis \
        --nproc ${task.cpus} \
        -o ${meta.id}.expected_cis.tsv \
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
    cooltools eigs-cis \
        --n-eigs ${params.compartment_n_eigs} \
        --phasing-track '${track}' \
        --out-prefix ${meta.id}.compartments \
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
    cooltools insulation \
        --nproc ${task.cpus} \
        --threshold ${params.insulation_threshold} \
        -o ${meta.id}.insulation.tsv \
        '${uri}' \
        ${windows}

    for window in ${windows}; do
        insulation_to_TAD.py \
            ${meta.id}.insulation.tsv \
            --window \${window} \
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
    cooltools dots \
        --nproc ${task.cpus} \
        --max-loci-separation ${params.max_loop_distance} \
        -o ${meta.id}.dots.tsv \
        '${uri}' \
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

    publishDir "${params.outdir}/06_hichip/peaks", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(bam), path(bai), path(input_peak)

    output:
    tuple val(meta), path("${meta.id}.peaks.bed"), emit: peaks

    when:
    params.run_hichip_peaks.toString().toBoolean()

    script:
    if (meta.peak_provided) {
        return """
        if [[ '${input_peak}' == *.gz ]]; then
            gzip -cd ${input_peak} | awk 'BEGIN{OFS="\\t"} !/^#/ && NF >= 3 {print}' > ${meta.id}.peaks.bed
        else
            awk 'BEGIN{OFS="\\t"} !/^#/ && NF >= 3 {print}' ${input_peak} > ${meta.id}.peaks.bed
        fi
        test -s ${meta.id}.peaks.bed
        """
    }
    if (!params.macs2_genome_size) {
        error "HiChIP sample ${meta.id} has no peak_file; set --macs2_genome_size (for example hs, mm, or an effective genome size)."
    }
    """
    macs2 callpeak \
        --treatment ${bam} \
        --format BAMPE \
        --gsize ${params.macs2_genome_size} \
        --name ${meta.id} \
        --qvalue ${params.macs2_qvalue} \
        --keep-dup all \
        ${params.macs2_extra_args}

    cp ${meta.id}_peaks.narrowPeak ${meta.id}.peaks.bed
    """

    stub:
    """
    touch ${meta.id}.peaks.bed
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
