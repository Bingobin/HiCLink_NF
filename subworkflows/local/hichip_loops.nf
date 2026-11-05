include { LOOP_INPUT_QC; CLOOPS_PRE; CLOOPS_CALL; CLOOPS_DIFF; MUSTACHE_PREPARE; MUSTACHE_CALL; MUSTACHE_DIFF } from '../../modules/local/hichip_loops'
include { FITHICHIP_PREPARE; FITHICHIP_UNION_PEAKS; FITHICHIP_LOOP_CALL; FITHICHIP_LOOP_DIFF } from '../../modules/local/fithichip_loops'


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
        fithichip         : [image: p.fithichip_loop_image.toString(), platform: p.fithichip_loop_platform.toString(),
                            qvalue: p.fithichip_loop_qvalue as Double,
                            background_fdr: p.fithichip_loop_background_fdr as Double,
                            diff_fdr: p.fithichip_loop_diff_fdr as Double,
                            fold_change: p.fithichip_loop_fold_change as Integer,
                            bcv: p.fithichip_loop_bcv as Double,
                            interaction_type: p.fithichip_loop_interaction_type as Integer,
                            p2p_background: p.fithichip_loop_p2p_background as Integer,
                            bias_type: p.fithichip_loop_bias_type as Integer],
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
    fitInputs = channel.empty()
    fitCalls = channel.empty()
    fitDiff = channel.empty()

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

    if ('fithichip' in cfg.callers) {
        fitHelper = file("${projectDir}/bin/fithichip.py", checkIfExists: true)
        fitLauncher = file("${projectDir}/bin/FitHiChIP_Docker.sh", checkIfExists: true)
        FITHICHIP_PREPARE(validated.map { meta, _pairs, matrix, _report ->
            tuple(meta, matrix, file(meta.peak_file, checkIfExists: true))
        }, cfg, [fitHelper, helper])
        fitInputs = FITHICHIP_PREPARE.out.inputs

        comparedIds = channel.value([ids: []])
        pairedCalls = channel.empty()
        if (cfg.run_diff) {
            comparedIds = checkedComparisons.collect().map { rows ->
                [ids: rows.collectMany { comp -> [comp.control, comp.treatment] }.unique()]
            }
            indexedFit = fitInputs.map { meta, inputs -> tuple(meta.id, meta, inputs, inputs) }
            pairedFit = pairLoopComparisons(checkedComparisons, indexedFit)
                .map { comp, treatment, _t, control, _c -> tuple(comp, treatment, control) }
            FITHICHIP_UNION_PEAKS(pairedFit, cfg, fitHelper)
            pairedCalls = FITHICHIP_UNION_PEAKS.out.groups.flatMap { comp, treatment, control, peaks ->
                // Refit both samples against the same peak universe in each comparison.
                [tuple([id: comp.control, run_id: "${comp.id}.control", role: 'control', comparison: comp], control, peaks),
                 tuple([id: comp.treatment, run_id: "${comp.id}.treatment", role: 'treatment', comparison: comp], treatment, peaks)]
            }
        }
        individualCalls = fitInputs.combine(comparedIds)
            .filter { meta, _inputs, used -> !used.ids.contains(meta.id) }
            .map { meta, inputs, _used ->
                tuple([id: meta.id, run_id: meta.id, role: 'individual'], inputs, inputs.resolve('peaks.bed'))
            }
        FITHICHIP_LOOP_CALL(pairedCalls.mix(individualCalls), cfg, fitHelper, fitLauncher)
        fitCalls = FITHICHIP_LOOP_CALL.out.calls

        if (cfg.run_diff) {
            fitDiffInput = fitCalls.filter { meta, _dir -> meta.role != 'individual' }
                .map { meta, dir -> tuple(meta.comparison.id, [meta: meta, dir: dir]) }
                .groupTuple(size: 2)
                .map { _id, rows ->
                    def control = rows.find { row -> row.meta.role == 'control' }
                    def treatment = rows.find { row -> row.meta.role == 'treatment' }
                    def comp = control.meta.comparison
                    tuple(comp, control.dir.resolve("${comp.control}.interactions_FitHiC.bed"),
                          treatment.dir.resolve("${comp.treatment}.interactions_FitHiC.bed"),
                          control.dir.resolve('chrom.sizes'))
                }
            FITHICHIP_LOOP_DIFF(fitDiffInput, cfg, fitHelper)
            fitDiff = FITHICHIP_LOOP_DIFF.out.differences
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
    fithichip_inputs = fitInputs
    fithichip_calls = fitCalls
    fithichip_diff = fitDiff
}
