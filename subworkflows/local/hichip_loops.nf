include { LOOP_INPUT_QC; CLOOPS_PRE; CLOOPS_CALL; CLOOPS_DIFF; MUSTACHE_PREPARE; MUSTACHE_CALL; MUSTACHE_DIFF } from '../../modules/local/hichip_loops'


def hichipLoopSettings(p) {
    return [
        callers           : p.hichip_loop_callers.toString().split(',').collect { name -> name.trim() }.findAll { name -> name }.unique(),
        run_diff          : p.hichip_loop_diff.toString().toBoolean(),
        resolution        : p.loop_resolution as Integer,
        min_distance      : p.min_loop_distance as Integer,
        max_distance      : p.max_loop_distance as Integer,
        chroms            : p.hichip_loop_chroms.toString().split(',').collect { chrom -> chrom.trim() },
        eps               : p.cloops2_eps.toString(),
        minpts            : p.cloops2_minpts.toString(),
        diff_pcut         : p.cloops2_diff_pcut as Double,
        diff_ma_fdr       : p.cloops2_diff_ma_fdr as Double,
        mustache_pt       : p.mustache_pt as Double,
        mustache_diff_pt  : p.mustache_diff_pt as Double,
        sparsity          : p.mustache_sparsity as Double,
        balance_max_iters : p.loop_balance_max_iters as Integer,
        outdir            : p.hichip_loop_outdir.toString()
    ]
}


def checkComparison(comparison, treatment, control) {
    if (comparison.treatment != treatment.id || comparison.control != control.id || treatment.id == control.id) {
        error "Invalid directed loop comparison: ${comparison}"
    }
    def tparts = treatment.id.tokenize('.')
    def cparts = control.id.tokenize('.')
    if (tparts.size() >= 3 && cparts.size() >= 3 && tparts.take(2) != cparts.take(2)) {
        error "Cross-cell-line or cross-IP loop comparison: ${treatment.id} vs ${control.id}"
    }
    return comparison
}


def pairLoopComparisons(comparisons, indexed) {
    // Both caller branches use [sample_id, meta, artifact1, artifact2].
    // combine supports repeated control keys; join would not provide a cross-product.
    return comparisons.map { comp -> tuple(comp.treatment, comp) }
        .combine(indexed, by: 0)
        .map { _id, comp, meta, first, second -> tuple(comp.control, comp, meta, first, second) }
        .combine(indexed, by: 0)
        .map { _id, comp, tmeta, tfirst, tsecond, cmeta, cfirst, csecond ->
            tuple(checkComparison(comp, tmeta, cmeta), tfirst, tsecond, cfirst, csecond)
        }
}


workflow HICHIP_LOOPS {
    take:
    samples          // queue: tuple(meta, pairs, matrix); meta.id and meta.assay are required
    comparisons      // queue: [id: comparison_id, treatment: sample_id, control: sample_id]
    cfg              // value: immutable settings map; no CSV parsing or local paths in this workflow

    main:
    helper = file("${projectDir}/bin/hichip_loops.py", checkIfExists: true)
    sampleList = samples.collect(flat: false).map { rows ->
        def ids = rows.collect { row -> row[0].id }
        if (!ids || ids.unique(false).size() != ids.size()) {
            error 'The HiChIP loop sample channel must contain nonempty, unique sample IDs.'
        }
        return rows
    }
    inputSamples = sampleList.flatMap { rows -> rows }
    sampleIndex = sampleList.map { rows -> [byId: rows.collectEntries { row -> [(row[0].id): row[0]] }] }
    checkedComparisons = comparisons.combine(sampleIndex).map { comp, index ->
        if (!(comp.id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) ||
            !index.byId.containsKey(comp.treatment) || !index.byId.containsKey(comp.control)) {
            error "Invalid comparison ID or missing comparison sample: ${comp}"
        }
        return checkComparison(comp, index.byId[comp.treatment], index.byId[comp.control])
    }
    LOOP_INPUT_QC(inputSamples, cfg, helper)
    validated = inputSamples
        .map { meta, pairs, matrix -> tuple(meta.id, meta, pairs, matrix) }
        .join(LOOP_INPUT_QC.out.reports.map { meta, report -> tuple(meta.id, report) })
        .map { _id, meta, pairs, matrix, report -> tuple(meta, pairs, matrix, report) }

    cloopsPre = channel.empty()
    cloopsCalls = channel.empty()
    cloopsDiff = channel.empty()
    mustacheCools = channel.empty()
    mustacheCalls = channel.empty()
    mustacheDiff = channel.empty()

    if ('cloops2' in cfg.callers) {
        CLOOPS_PRE(validated.map { meta, pairs, _matrix, report -> tuple(meta, pairs, report) }, cfg, helper)
        CLOOPS_CALL(CLOOPS_PRE.out.pre, cfg, helper)
        cloopsPre = CLOOPS_PRE.out.pre
        cloopsCalls = CLOOPS_CALL.out.calls

        if (cfg.run_diff) {
            indexed = CLOOPS_PRE.out.pre.map { meta, pre -> tuple(meta.id, meta, pre) }
                .join(CLOOPS_CALL.out.calls.map { meta, calls -> tuple(meta.id, calls.resolve("${meta.id}_loops.txt")) })
            diffInput = pairLoopComparisons(checkedComparisons, indexed)
            CLOOPS_DIFF(diffInput, cfg, helper)
            cloopsDiff = CLOOPS_DIFF.out.differences
        }
    }

    if ('mustache' in cfg.callers) {
        MUSTACHE_PREPARE(validated.map { meta, _pairs, matrix, report -> tuple(meta, matrix, report) }, cfg, helper)
        MUSTACHE_CALL(MUSTACHE_PREPARE.out.cools, cfg, helper)
        mustacheCools = MUSTACHE_PREPARE.out.cools
        mustacheCalls = MUSTACHE_CALL.out.calls

        if (cfg.run_diff) {
            indexedCools = MUSTACHE_PREPARE.out.cools.map { meta, cool, report -> tuple(meta.id, meta, cool, report) }
            diffCoolInput = pairLoopComparisons(checkedComparisons, indexedCools)
            MUSTACHE_DIFF(diffCoolInput, cfg, helper)
            mustacheDiff = MUSTACHE_DIFF.out.differences
        }
    }

    emit:
    input_qc = LOOP_INPUT_QC.out.reports
    cloops_pre = cloopsPre
    cloops_calls = cloopsCalls
    cloops_diff = cloopsDiff
    mustache_cools = mustacheCools
    mustache_calls = mustacheCalls
    mustache_diff = mustacheDiff
}
