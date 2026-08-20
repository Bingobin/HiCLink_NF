#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "$projectDir/assets/samplesheet.csv"
params.outdir = "results"
params.ref = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta"
params.chrom = "/lustre/home/acct-medty/medty-c/script/wgs_pipeline/hg38.chrom.size"


log.info """\
    H I C    --    P I P E L I N E
    ===================================
    Sample Info :  ${params.input}
    Out Dir     :  ${params.outdir}
    Reference   :  ${params.ref}
    Chrom.size  :  ${params.chrom}
    """
    .stripIndent(true)

process BWA_ALIGN {
    tag "bwa alignment in $ID"

    input:
    tuple val(ID), val(R1), val(R2)

    output:
    tuple val(ID), path("${ID}.aln.bam")

    script:
    """
    bwa mem -SP -t $task.cpus ${params.ref}  $R1 $R2 | samtools view  -@ $task.cpus -bS - > ${ID}.aln.bam
    """
}

process PAIRTOOLS {
    tag "pairtools in $ID"
    publishDir "$params.outdir", pattern: "${ID}.*", mode: 'copy'

    input:
    tuple val(ID), path(BAM)

    output:
    tuple val(ID), path("${ID}.dedup.pairs"), emit: pairs
    tuple val(ID), path("${ID}.*.stats"), emit: stats

    script:
    """
    pairtools parse -c ${params.chrom} -o - --drop-sam $BAM | pairtools sort -o - | pairtools dedup -o ${ID}.dedup.pairs
    pairtools stats -o ${ID}.dedup.pairs.stats ${ID}.dedup.pairs
    """
}

process BAMstats {
    tag "samtools stats in $ID"
    publishDir "$params.outdir", pattern: "${ID}.aln.bam.stats", mode: 'copy'

    input:
    tuple val(ID), path(BAM)

    output:
    tuple val(ID), path("${ID}.*.stats")

    script:
    """
    samtools stats -@ $task.cpus ${BAM} > ${ID}.aln.bam.stats
    """
}

process COOLER {
    tag "cooler in $ID"
    publishDir "$params.outdir", pattern: "${ID}.*cool", mode: 'copy'


    input:
    tuple val(ID), path(PAIR)

    output:
    tuple val(ID), path("${ID}.10kb.cool"), emit: cool
    tuple val(ID), path("${ID}.10kb.mcool"), emit: mcool

    script:
    """
    cooler makebins ${params.chrom} 10000 > bins_10kb.bed
    cooler cload pairs --assembly hg38 -c1 2 -p1 3 -c2 4 -p2 5 bins_10kb.bed $PAIR ${ID}.10kb.cool
    cooler balance ${ID}.10kb.cool
    cooler zoomify -r N,250000 ${ID}.10kb.cool
    """
}

workflow{
    Channel.fromPath(params.input)
    .splitCsv(header:true,sep:',')
    .map { ["${it.ID}","${it.R1}","${it.R2}"] }
    .set{ch_fq}
    ch_aln = BWA_ALIGN(ch_fq)
    ch_pair = PAIRTOOLS(ch_aln)
//    BAMstats(ch_aln)
    ch_cool = COOLER(ch_pair.pairs)

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. someting went wrong" )
}
