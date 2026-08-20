#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "$projectDir/assets/samplesheet.cool.csv"
params.outdir = "results_proc"
params.ref = "/lustre/home/acct-medkkw/medlyb/database/annotation/gatk_ann/hg38/bwaindex2/Homo_sapiens_assembly38.fasta"
params.chrom = "/lustre/home/acct-medty/medty-c/script/wgs_pipeline/hg38.chrom.size"
params.res = 10000


log.info """\
    HIC process fo compare -- PIPELINE
    ===================================
    Sample Info :  ${params.input}
    Out Dir     :  ${params.outdir}
    Reference   :  ${params.ref}
    Chrom.size  :  ${params.chrom}
    Resolution  :  ${params.res}
    """
    .stripIndent(true)

process COOL_QC {
    cpus 2
    
    tag "COOL QC in $ID"
    publishDir "$params.outdir/QC", pattern: "${ID}*.txt", mode: 'copy'

    input:
    tuple val(ID), val(MCOOL), val(CTR_MCOOL)

    output:
    tuple val(ID), path("${ID}*.txt")

    script:
    """
    cooler balance ${MCOOL}::/resolutions/${params.res * 10} --force
    cooltools coverage ${MCOOL}::/resolutions/${params.res} -o  ${ID}.${params.res / 1000}kb.coverage.txt
    cooltools coverage ${MCOOL}::/resolutions/${params.res * 10} -o  ${ID}.${params.res / 100}kb.coverage.txt
    cool_sparsity.sh ${MCOOL}::/resolutions/${params.res} > ${ID}.${params.res / 1000}kb.sparsity.txt
    cool_sparsity.sh ${MCOOL}::/resolutions/${params.res * 10} > ${ID}.${params.res / 100}kb.sparsity.txt
    """
}

process AB_Compartments {
    cpus 2

    tag "A/B compart in $ID"
    publishDir "$params.outdir/AB_Compartments", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), val(MCOOL), val(CTR_MCOOL)

    output:
    tuple val(ID), path("${ID}*")

    script:
    """
    cooltools eigs-cis ${MCOOL}::/resolutions/${params.res * 10} --phasing-track $projectDir/assets/HG38_GC_${params.res / 100}kb.bedgraph --n-eigs 3 --out-prefix ${ID}.${params.res / 100}kb.eigs
    """
}

process COOL_TAD {
    cpus 2

    tag "cooltools TAD in $ID"
    publishDir "$params.outdir/TAD", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), val(MCOOL), val(CTR_MCOOL)

    output:
    tuple val(ID), path("${ID}*")

    script:
    """
    cooltools insulation ${MCOOL}::/resolutions/${params.res} ${params.res * 5} ${params.res * 10} ${params.res * 25} ${params.res * 50} > ${ID}.${params.res / 1000}kb.insulation.tsv
    insulation_to_TAD.py -o ${ID}.${params.res * 10 / 1000}kb.TAD.bed ${ID}.${params.res / 1000}kb.insulation.tsv -w ${params.res * 10}
    insulation_to_TAD.py -o ${ID}.${params.res * 25 / 1000}kb.TAD.bed ${ID}.${params.res / 1000}kb.insulation.tsv -w ${params.res * 25}
    """
}

process DIFF_LOOP {
    cpus 2

    tag "diff_mustache.py in $ID"
    publishDir "$params.outdir/LOOP", pattern: "${ID}*", mode: 'copy'

    input:
    tuple val(ID), val(MCOOL), val(CTR_MCOOL)

    output:
    tuple val(ID), path("${ID}*")

    script:
    """
    cooler cp ${MCOOL}::/resolutions/${params.res} ${ID}_t.${params.res / 1000}kb.cool
    cooler cp ${CTR_MCOOL}::/resolutions/${params.res} ${ID}_c.${params.res / 1000}kb.cool
    python3 ~/soft/mustache/mustache/diff_mustache.py -f1 ${ID}_t.${params.res / 1000}kb.cool -f2 ${ID}_c.${params.res / 1000}kb.cool -pt 0.05 -pt2 0.1 -o ${ID}.${params.res / 1000}kb -r ${params.res} -st 0.8 1>${ID}.${params.res / 1000}kb.mustache.out 2>${ID}.${params.res / 1000}kb.mustache.err
    awk 'BEGIN{getline}{OFS="\\t";print \$1,\$2,\$3,\$4,\$5,\$6,".",\$7,"+","+"}' ${ID}.${params.res / 1000}kb.diffloop1  > ${ID}.${params.res / 1000}kb.diffloop1.bedpe
    awk 'BEGIN{getline}{OFS="\\t";print \$1,\$2,\$3,\$4,\$5,\$6,".",\$7,"+","+"}' ${ID}.${params.res / 1000}kb.diffloop2  > ${ID}.${params.res / 1000}kb.diffloop2.bedpe
    """
}

workflow{
    Channel.fromPath(params.input)
    .splitCsv(header:true,sep:',')
    .map { ["${it.ID}","${it.MCOOL}","${it.CTR_MCOOL}"] }
    .set{ch_cool}
    COOL_QC(ch_cool)
    AB_Compartments(ch_cool)
    COOL_TAD(ch_cool)

    ch_cool.filter { it[2] != 'NA'}.set{ch_paired_cool}
//    ch_paired_cool.view()
    DIFF_LOOP(ch_paired_cool)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. someting went wrong" )
}
