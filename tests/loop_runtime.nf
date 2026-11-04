#!/usr/bin/env nextflow
// Harmless resource/configuration regression: neither process executes a loop caller.
nextflow.enable.dsl = 2

process CLOOPS_CALL {
    label 'hichip_loop_analysis'
    output:
    path 'cloops_cpus.txt'
    script:
    """
    printf '%s\\n' '${task.cpus}' > cloops_cpus.txt
    """
}

process MUSTACHE_CALL {
    label 'hichip_loop_analysis'
    output:
    path 'mustache_cpus.txt'
    script:
    """
    printf '%s\\n' '${task.cpus}' > mustache_cpus.txt
    command -v python >> mustache_cpus.txt
    printf '%s\\n' '${task.memory.toBytes()}' >> mustache_cpus.txt
    """
}

workflow {
    CLOOPS_CALL()
    MUSTACHE_CALL()
    CLOOPS_CALL.out.combine(MUSTACHE_CALL.out).view { cloops, mustache ->
        def expected = System.getProperty('os.name').toLowerCase().contains('mac')
            ? 1 : (params.hichip_loop_cpus as Integer)
        def actual = mustache.readLines()
        if (cloops.text.trim() != expected.toString() ||
            actual[0] != params.hichip_loop_cpus.toString()) {
            error "CPU routing failed: cLoops2=${cloops.text}, Mustache=${mustache.text}"
        }
        if (params.hichip_loops_env && actual[1] != "${params.hichip_loops_env}/bin/python") {
            error "Local environment prefix was not applied: ${actual[1]}"
        }
        def expectedMemory = (params.memory_per_cpu_gb as Long) * (params.hichip_loop_cpus as Long) * 1024L * 1024L * 1024L
        if ((actual[2] as Long) != expectedMemory) {
            error "Late-bound memory settings were not applied: ${actual[2]} vs ${expectedMemory}"
        }
        return "Runtime regression passed: cLoops2=${expected}, Mustache=${params.hichip_loop_cpus}; no analysis tools executed."
    }
}
