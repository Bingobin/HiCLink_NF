process VALIDATE_SAMPLESHEET {
    tag 'samplesheet'
    label 'small'

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode, overwrite: true

    input:
    path input_sheet
    val base_dir

    output:
    path 'samplesheet.validated.csv', emit: csv

    script:
    """
    validate_samplesheet.py \\
        --input ${input_sheet} \\
        --output samplesheet.validated.csv \\
        --base-dir '${base_dir}'
    """

    stub:
    """
    cp ${input_sheet} samplesheet.validated.csv
    """
}


process FASTQC_RAW {
    tag "${meta.id}"
    label 'small'

    publishDir "${params.outdir}/01_qc/fastqc_raw", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}.*_fastqc.zip"), emit: zips
    tuple val(meta), path("${meta.id}.*_fastqc.html"), emit: html

    when:
    params.run_fastqc.toString().toBoolean()

    script:
    """
    fastqc --threads ${task.cpus} ${read1} ${read2}
    for report in *_fastqc.zip *_fastqc.html; do
        mv "\${report}" "${meta.id}.\${report}"
    done
    """

    stub:
    """
    touch ${meta.id}.${read1.simpleName}_fastqc.zip ${meta.id}.${read2.simpleName}_fastqc.zip
    touch ${meta.id}.${read1.simpleName}_fastqc.html ${meta.id}.${read2.simpleName}_fastqc.html
    """
}


process FASTP {
    tag "${meta.id}"
    label 'large'

    publishDir "${params.outdir}/01_qc/fastp", pattern: '*.{json,html}', mode: params.publish_mode, overwrite: true
    publishDir "${params.outdir}/02_preprocessing/trimmed_fastq", pattern: '*.fastq.gz', mode: params.publish_mode, overwrite: true, enabled: params.save_trimmed.toString().toBoolean()

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.id}.R1.fastq.gz"), path("${meta.id}.R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json
    tuple val(meta), path("${meta.id}.fastp.html"), emit: html

    script:
    def noTrim = params.skip_trimming.toString().toBoolean()
        ? '--disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g'
        : params.fastp_extra_args
    """
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${meta.id}.R1.fastq.gz \\
        --out2 ${meta.id}.R2.fastq.gz \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        ${noTrim}
    """

    stub:
    """
    touch ${meta.id}.R1.fastq.gz ${meta.id}.R2.fastq.gz
    echo '{}' > ${meta.id}.fastp.json
    touch ${meta.id}.fastp.html
    """
}


process BWA_ALIGN {
    tag "${meta.id}"
    label 'large'

    publishDir "${params.outdir}/01_qc/alignment", pattern: '*.flagstat', mode: params.publish_mode, overwrite: true
    publishDir "${params.outdir}/02_preprocessing/alignment_bam", pattern: '*.{peaks.bam,peaks.bam.bai}', mode: params.publish_mode, overwrite: true, enabled: params.save_alignment_bam.toString().toBoolean()

    input:
    tuple val(meta), path(read1), path(read2), path(fasta), path(bwa_indexes)

    output:
    tuple val(meta), path("${meta.id}.contacts.bam"), emit: contacts
    tuple val(meta), path("${meta.id}.peaks.bam"), path("${meta.id}.peaks.bam.bai"), emit: peak_bam
    tuple val(meta), path("${meta.id}.contacts.flagstat"), path("${meta.id}.peaks.flagstat"), emit: flagstats

    script:
    """
    set -o pipefail

    bwa mem -SP -t ${task.cpus} ${fasta} ${read1} ${read2} \\
        | samtools view -@ ${task.cpus} -bh -o ${meta.id}.contacts.bam -

    samtools flagstat -@ ${task.cpus} ${meta.id}.contacts.bam > ${meta.id}.contacts.flagstat

    samtools sort -n -@ ${task.cpus} -o ${meta.id}.namesort.bam ${meta.id}.contacts.bam
    samtools fixmate -m -@ ${task.cpus} ${meta.id}.namesort.bam ${meta.id}.fixmate.bam
    rm -f ${meta.id}.namesort.bam

    samtools sort -@ ${task.cpus} -o ${meta.id}.positionsort.bam ${meta.id}.fixmate.bam
    rm -f ${meta.id}.fixmate.bam

    samtools markdup -r -@ ${task.cpus} ${meta.id}.positionsort.bam ${meta.id}.markdup.bam
    rm -f ${meta.id}.positionsort.bam

    samtools view -@ ${task.cpus} -bh -q ${params.min_mapq} -F 3844 \\
        -o ${meta.id}.peaks.bam ${meta.id}.markdup.bam
    rm -f ${meta.id}.markdup.bam

    samtools index -@ ${task.cpus} ${meta.id}.peaks.bam
    samtools flagstat -@ ${task.cpus} ${meta.id}.peaks.bam > ${meta.id}.peaks.flagstat
    """

    stub:
    """
    touch ${meta.id}.contacts.bam
    touch ${meta.id}.peaks.bam ${meta.id}.peaks.bam.bai
    touch ${meta.id}.contacts.flagstat ${meta.id}.peaks.flagstat
    """
}


process PEAKS_BIGWIG {
    tag "${meta.id}"
    label 'medium'

    publishDir "${params.outdir}/06_hichip/bigwig", pattern: '*.peaks.RPGC.bw', mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.peaks.RPGC.bw"), emit: bigwig

    when:
    params.save_peak_bigwig.toString().toBoolean() && meta.assay == 'hichip'

    script:
    def genomeKey = params.genome?.toString()?.trim()?.toLowerCase()
    def defaultEffectiveGenomeSizes = [
        'hg38'  : 2913022398L,
        'grch38': 2913022398L
    ]
    def effectiveGenomeSizeValue = params.effective_genome_size ?: defaultEffectiveGenomeSizes[genomeKey]
    if (!effectiveGenomeSizeValue) {
        error "HiChIP BigWig generation requires --effective_genome_size for genome '${params.genome}'. hg38/GRCh38 defaults to 2913022398; custom genomes must set this parameter explicitly."
    }
    long effectiveGenomeSize
    int binSize
    try {
        effectiveGenomeSize = effectiveGenomeSizeValue as Long
        binSize = params.bigwig_bin_size as Integer
    } catch (Exception ignored) {
        error "Parameters --effective_genome_size and --bigwig_bin_size must be positive integers."
    }
    if (effectiveGenomeSize <= 0 || binSize <= 0) {
        error "Parameters --effective_genome_size and --bigwig_bin_size must be positive integers."
    }
    def ignoredChromosomes = params.bigwig_ignore_for_normalization?.toString()?.trim()
    def ignoreArgs = ignoredChromosomes
        ? "--ignoreForNormalization ${ignoredChromosomes.split(',').collect { it.trim() }.findAll { it }.join(' ')}"
        : ''
    """
    # Do not use --extendReads: HiChIP mates represent proximity contacts, not
    # contiguous ChIP-seq fragments.
    bamCoverage \\
        --bam ${bam} \\
        --outFileName ${meta.id}.peaks.RPGC.bw \\
        --outFileFormat bigwig \\
        --binSize ${binSize} \\
        --normalizeUsing RPGC \\
        --effectiveGenomeSize ${effectiveGenomeSize} \\
        --numberOfProcessors ${task.cpus} \\
        ${ignoreArgs} \\
        ${params.bigwig_extra_args}
    """

    stub:
    """
    touch ${meta.id}.peaks.RPGC.bw
    """
}


process PAIRTOOLS {
    tag "${meta.id}"
    label 'medium'

    publishDir "${params.outdir}/03_pairs", pattern: '*.valid.pairs.gz', mode: params.publish_mode, overwrite: true, enabled: params.save_pairs.toString().toBoolean()
    publishDir "${params.outdir}/01_qc/pairtools", pattern: '*.stats', mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(bam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.id}.valid.pairs.gz"), emit: pairs
    tuple val(meta), path("${meta.id}.parse.stats"), path("${meta.id}.dedup.stats"), path("${meta.id}.pairtools.stats"), emit: stats

    script:
    """
    set -o pipefail

    pairtools parse \\
        --chroms-path ${chrom_sizes} \\
        --min-mapq ${params.min_mapq} \\
        --walks-policy ${params.walks_policy} \\
        --drop-sam \\
        --output-stats ${meta.id}.parse.stats \\
        --nproc-in ${task.cpus} \\
        --nproc-out ${task.cpus} \\
        ${bam} \\
    | pairtools sort \\
        --nproc ${task.cpus} \\
        --tmpdir . \\
    | pairtools dedup \\
        --max-mismatch ${params.dedup_max_mismatch} \\
        --output-stats ${meta.id}.dedup.stats \\
    | pairtools select \\
        '((pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU"))' \\
        -o ${meta.id}.valid.pairs.gz

    pairtools stats -o ${meta.id}.pairtools.stats ${meta.id}.valid.pairs.gz
    """

    stub:
    """
    touch ${meta.id}.valid.pairs.gz
    touch ${meta.id}.parse.stats ${meta.id}.dedup.stats ${meta.id}.pairtools.stats
    """
}


process COOLER_MATRIX {
    tag "${meta.id}"
    label 'large'

    publishDir "${params.outdir}/04_matrices", mode: params.publish_mode, overwrite: true

    input:
    tuple val(meta), path(pairs), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.id}.${params.base_resolution}.cool"), path("${meta.id}.mcool"), emit: matrices

    script:
    """
    cooler cload pairs \\
        --assembly ${params.genome} \\
        -c1 2 -p1 3 -c2 4 -p2 5 \\
        ${chrom_sizes}:${params.base_resolution} \\
        ${pairs} \\
        ${meta.id}.${params.base_resolution}.cool

    cooler balance \\
        --nproc ${task.cpus} \\
        --ignore-diags ${params.balance_ignore_diags} \\
        ${meta.id}.${params.base_resolution}.cool

    cooler zoomify \\
        --balance \\
        --balance-args '--nproc ${task.cpus} --ignore-diags ${params.balance_ignore_diags}' \\
        --nproc ${task.cpus} \\
        --resolutions ${params.resolutions} \\
        -o ${meta.id}.mcool \\
        ${meta.id}.${params.base_resolution}.cool
    """

    stub:
    """
    touch ${meta.id}.${params.base_resolution}.cool ${meta.id}.mcool
    """
}
