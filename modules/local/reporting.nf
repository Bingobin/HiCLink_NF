process MULTIQC {
    tag 'project'
    label 'small'

    publishDir "${params.outdir}/pipeline_info/multiqc", mode: params.publish_mode, overwrite: true

    input:
    path reports
    path multiqc_config

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data

    when:
    params.run_multiqc.toString().toBoolean()

    script:
    """
    multiqc \
        --force \
        --config ${multiqc_config} \
        --filename multiqc_report.html \
        --outdir . \
        .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html multiqc_data/multiqc_stub.txt
    """
}
