def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}


process DISCOVER_LOOP_INPUTS {
    tag 'published-input-discovery'
    label 'hichip_loop_qc'
    // Recheck lightweight discovery/signatures even when analysis uses -resume.
    cache false
    publishDir { "${cfg.outdir}/pipeline_info" }, mode: 'copy', overwrite: false

    input:
    path sheet
    val published_root
    val cfg
    path helper

    output:
    path 'loop_samples.tsv', emit: samples
    path 'loop_comparisons.tsv', emit: comparisons
    path 'loop_settings.json', emit: settings
    path 'loop_inputs_signature.json', emit: signature

    script:
    def settings = shellQuote(groovy.json.JsonOutput.toJson(cfg))
    """
    python '${helper}' discover --samplesheet '${sheet}' \\
        --results-dir ${shellQuote(published_root)} --settings ${settings} --output-dir .
    """
}


process LOOP_INPUT_QC {
    tag "${meta.id}"
    label 'hichip_loop_qc'
    publishDir { "${cfg.outdir}/input_qc" }, mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(pairs), path(matrix)
    val cfg
    path helper

    output:
    tuple val(meta), path("${meta.id}.input_qc.json"), emit: reports

    script:
    """
    python '${helper}' check-sample --meta ${shellQuote(groovy.json.JsonOutput.toJson(meta))} \\
        --pairs '${pairs}' --matrix '${matrix}' --settings ${shellQuote(groovy.json.JsonOutput.toJson(cfg))} \\
        --output '${meta.id}.input_qc.json'
    """
}


process CLOOPS_PRE {
    tag "${meta.id}"
    label 'hichip_loop_analysis'
    publishDir { "${cfg.outdir}/cloops2/pre" }, mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(pairs), path(qc)
    val cfg
    path helper

    output:
    tuple val(meta), path("${meta.id}"), emit: pre

    script:
    """
    cLoops2 pre -f '${pairs}' -format pairs -o '${meta.id}' \\
        -c '${cfg.chroms.join(',')}' -cut ${cfg.min_distance} -mcut ${cfg.max_distance} -p ${task.cpus}
    python '${helper}' check-pre '${meta.id}/petMeta.json'
    """
}


process CLOOPS_CALL {
    tag "${meta.id}"
    label 'hichip_loop_analysis'
    // cLoops2 metadata contains absolute paths: never update a staged symlink.
    stageInMode 'copy'
    publishDir { "${cfg.outdir}/cloops2/call" }, mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(pre_dir, stageAs: 'pre_input')
    val cfg
    path helper

    output:
    tuple val(meta), path("${meta.id}"), emit: calls

    script:
    """
    cLoops2 update -d '${pre_dir}'
    python '${helper}' check-pre '${pre_dir}/petMeta.json'
    mkdir '${meta.id}'
    cLoops2 callLoops -d '${pre_dir}' -o '${meta.id}/${meta.id}' \\
        -eps '${cfg.eps}' -minPts '${cfg.minpts}' -hic \\
        -cut ${cfg.min_distance} -mcut ${cfg.max_distance} -p ${task.cpus} -w -j
    test -s '${meta.id}/${meta.id}_loops.txt'
    """
}


process CLOOPS_DIFF {
    tag "${comparison.id}"
    label 'hichip_loop_analysis'
    stageInMode 'copy'
    publishDir { "${cfg.outdir}/cloops2/diff" }, mode: 'copy', overwrite: false

    input:
    tuple val(comparison), path(tpre, stageAs: 'treatment_pre'), path(tloop, stageAs: 'treatment_loops.txt'), path(cpre, stageAs: 'control_pre'), path(cloop, stageAs: 'control_loops.txt')
    val cfg
    path helper

    output:
    tuple val(comparison), path("${comparison.id}"), emit: differences

    script:
    """
    cLoops2 update -d '${tpre}'
    cLoops2 update -d '${cpre}'
    python '${helper}' require-loops '${tloop}' '${cloop}'
    mkdir '${comparison.id}'
    printf 'role\\tsample_id\\ntreatment\\t%s\\ncontrol\\t%s\\n' \\
        '${comparison.treatment}' '${comparison.control}' > '${comparison.id}/comparison.tsv'
    cLoops2 callDiffLoops -tloop '${tloop}' -cloop '${cloop}' -td '${tpre}' -cd '${cpre}' \\
        -o '${comparison.id}/${comparison.id}' -cut ${cfg.min_distance} -mcut ${cfg.max_distance} \\
        -pcut ${cfg.diff_pcut} -fdr ${cfg.diff_ma_fdr} -p ${task.cpus} -w -j
    test -s '${comparison.id}/${comparison.id}_dloops.txt'
    """
}


process MUSTACHE_PREPARE {
    tag "${meta.id}"
    label 'hichip_loop_analysis'
    publishDir { "${cfg.outdir}/mustache/cool" }, pattern: '*.cool', mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(matrix), path(qc)
    val cfg
    path helper

    output:
    tuple val(meta), path("${meta.id}.${cfg.resolution}.cool"), path(qc), emit: cools

    script:
    """
    python '${helper}' prepare-cool --source '${matrix}' --output '${meta.id}.${cfg.resolution}.cool' \\
        --settings ${shellQuote(groovy.json.JsonOutput.toJson(cfg))} --cpus ${task.cpus}
    """
}


process MUSTACHE_CALL {
    tag "${meta.id}"
    label 'hichip_loop_analysis'
    publishDir { "${cfg.outdir}/mustache/call" }, mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(cool), path(qc)
    val cfg
    path helper

    output:
    tuple val(meta), path("${meta.id}"), emit: calls

    script:
    def chroms = cfg.chroms.collect { chrom -> shellQuote(chrom) }.join(' ')
    """
    mkdir '${meta.id}'
    mustache -f '${cool}' -r ${cfg.resolution} -d ${cfg.max_distance} \\
        -pt ${cfg.mustache_pt} -st ${cfg.sparsity} -p ${task.cpus} -norm weight \\
        -ch ${chroms} -o '${meta.id}/${meta.id}.native.tsv'
    python '${helper}' filter-mustache '${meta.id}/${meta.id}.native.tsv' \\
        '${meta.id}/${meta.id}.loops.tsv' ${cfg.min_distance} ${cfg.max_distance}
    """
}


process MUSTACHE_DIFF {
    tag "${comparison.id}"
    label 'hichip_loop_analysis'
    publishDir { "${cfg.outdir}/mustache/diff" }, mode: 'copy', overwrite: false

    input:
    tuple val(comparison), path(tcool, stageAs: 'treatment.cool'), path(tqc, stageAs: 'treatment_qc.json'), path(ccool, stageAs: 'control.cool'), path(cqc, stageAs: 'control_qc.json')
    val cfg
    path helper

    output:
    tuple val(comparison), path("${comparison.id}"), emit: differences

    script:
    def chroms = cfg.chroms.collect { chrom -> shellQuote(chrom) }.join(' ')
    """
    python '${helper}' check-pair-reports '${tqc}' '${cqc}'
    mkdir '${comparison.id}'
    # Mustache 1.3.3 requires direct execution of the package file, not python -m.
    diff_script=\$(python -c 'import importlib.util,pathlib; print(pathlib.Path(importlib.util.find_spec("mustache").origin).parent / "diff_mustache.py")')
    python "\$diff_script" -f1 '${tcool}' -f2 '${ccool}' -r ${cfg.resolution} -d ${cfg.max_distance} \\
        -pt ${cfg.mustache_pt} -pt2 ${cfg.mustache_diff_pt} -st ${cfg.sparsity} -p ${task.cpus} \\
        -norm weight -ch ${chroms} -o '${comparison.id}/${comparison.id}'
    for suffix in loop1 loop2 diffloop1 diffloop2; do
        case "\$suffix" in
            loop1) label=treatment.loops ;;
            loop2) label=control.loops ;;
            diffloop1) label=treatment_enriched.candidates ;;
            diffloop2) label=control_enriched.candidates ;;
        esac
        python '${helper}' filter-mustache '${comparison.id}/${comparison.id}.'"\$suffix" \\
            '${comparison.id}/${comparison.id}.'"\$label"'.tsv' ${cfg.min_distance} ${cfg.max_distance}
    done
    """
}
