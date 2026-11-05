def fitSettingsJson(cfg) {
    return "'" + groovy.json.JsonOutput.toJson(cfg).replace("'", "'\"'\"'") + "'"
}

process FITHICHIP_PREPARE {
    tag "${meta.id}"
    label 'hichip_loop_analysis'
    publishDir { "${cfg.outdir}/fithichip/inputs" }, mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(matrix), path(peaks)
    val cfg
    path helpers

    output:
    tuple val(meta), path("${meta.id}.fit_inputs"), emit: inputs

    script:
    """
    python fithichip.py prepare --matrix '${matrix}' --peaks '${peaks}' \\
        --output '${meta.id}.fit_inputs' --cpus ${task.cpus} \\
        --settings ${fitSettingsJson(cfg)}
    """
}

process FITHICHIP_UNION_PEAKS {
    tag "${comparison.id}"
    label 'hichip_loop_qc'
    publishDir { "${cfg.outdir}/fithichip/peaks/${comparison.id}" }, pattern: '*.bed', mode: 'copy', overwrite: false

    input:
    tuple val(comparison), path(treatment_inputs), path(control_inputs)
    val cfg
    path helper

    output:
    tuple val(comparison), path(treatment_inputs), path(control_inputs), path('union.peaks.bed'), emit: groups

    script:
    """
    cmp '${control_inputs}/chrom.sizes' '${treatment_inputs}/chrom.sizes'
    python '${helper}' union '${control_inputs}/peaks.bed' '${treatment_inputs}/peaks.bed' \\
        --output union.peaks.bed
    """
}

process FITHICHIP_LOOP_CALL {
    tag "${meta.run_id}"
    label 'hichip_loop_fithichip'
    publishDir { "${cfg.outdir}/fithichip/call" }, mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(inputs), path(peaks, stageAs: 'anchors/peaks.bed')
    val cfg
    path helper
    path launcher

    output:
    tuple val(meta), path("${meta.run_id}.fithichip"), emit: calls

    script:
    """
    python '${helper}' config --inputs '${inputs}' --peaks '${peaks}' --sample '${meta.id}' \\
        --output-dir '${meta.run_id}.fithichip' --output sample.config.txt \\
        --settings ${fitSettingsJson(cfg)}
    bash '${launcher}' -C sample.config.txt --image '${cfg.fithichip.image}' \\
        --platform '${cfg.fithichip.platform}' --cpus ${task.cpus} --memory-bytes ${task.memory.toBytes()}
    """
}

process FITHICHIP_LOOP_DIFF {
    tag "${comparison.id}"
    label 'hichip_loop_fithichip'
    publishDir { "${cfg.outdir}/fithichip/diff" }, mode: 'copy', overwrite: false

    input:
    tuple val(comparison), path(control_loops, stageAs: 'control.bed'), path(treatment_loops, stageAs: 'treatment.bed'), path(chrom_sizes, stageAs: 'chrom.sizes')
    val cfg
    path helper

    output:
    tuple val(comparison), path("${comparison.id}.fithichip_diff"), emit: differences

    script:
    """
    python '${helper}' diff --control '${control_loops}' --treatment '${treatment_loops}' \\
        --chrom-sizes '${chrom_sizes}' --output '${comparison.id}.fithichip_diff' \\
        --comparison '${groovy.json.JsonOutput.toJson(comparison)}' \\
        --settings ${fitSettingsJson(cfg)} --image '${cfg.fithichip.image}' \\
        --platform '${cfg.fithichip.platform}' --cpus ${task.cpus} --memory-bytes ${task.memory.toBytes()}
    """
}
