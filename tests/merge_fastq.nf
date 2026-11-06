#!/usr/bin/env nextflow
// Execute the real merge process on tiny synthetic clean FASTQs, without other tools.
nextflow.enable.dsl = 2
include { MERGE_FASTQ } from '../modules/local/preprocessing'

workflow {
    inputs = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.ID, row.run_id, file(row.R1), file(row.R2)) }
        .groupTuple()
        .map { id, runs, r1, r2 -> tuple([id: id], runs, r1, r2) }
    MERGE_FASTQ(inputs)
    MERGE_FASTQ.out.reads.view { meta, r1, r2 -> "MERGED\t${meta.id}\t${r1}\t${r2}" }
}
