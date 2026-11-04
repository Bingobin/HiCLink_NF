#!/usr/bin/env nextflow
// Routing-only test: no process, tool, or sequencing data is executed/read.
nextflow.enable.dsl = 2
include { pairLoopComparisons } from '../subworkflows/local/hichip_loops'

workflow {
    indexed = channel.of('CTRL', 'SH1', 'SH2')
        .map { id -> tuple(id, [id: id], "${id}.artifact1", "${id}.artifact2") }
    comparisons = channel.of(
        [id: 'SH1_vs_CTRL', treatment: 'SH1', control: 'CTRL'],
        [id: 'SH2_vs_CTRL', treatment: 'SH2', control: 'CTRL']
    )
    pairLoopComparisons(comparisons, indexed).collect(flat: false).view { rows ->
        if (rows.size() != 2 || rows.collect { row -> row[0].id }.sort() != ['SH1_vs_CTRL', 'SH2_vs_CTRL']) {
            error "Incorrect shared-control routing: ${rows}"
        }
        rows.each { row ->
            if (row[1] != "${row[0].treatment}.artifact1" || row[3] != 'CTRL.artifact1') {
                error "Incorrect treatment/control artifact assignment: ${row}"
            }
        }
        return 'Routing-only test passed: two treatments matched to the same control; no analysis tasks.'
    }
}
