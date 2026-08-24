process MULTIQC {
    tag 'project'
    label 'small'

    publishDir "${params.outdir}/pipeline_info/multiqc", mode: params.publish_mode, overwrite: true

    input:
    path reports
    path multiqc_config

    output:
    path 'multiqc_report.html', emit: report
    path '*_data', emit: data
    path '*_plots', optional: true, emit: plots

    when:
    params.run_multiqc.toString().toBoolean()

    script:
    """
    multiqc \\
        --force \\
        --config ${multiqc_config} \\
        --filename multiqc_report.html \\
        --outdir . \\
        .
    """

    stub:
    """
    mkdir -p multiqc_report_data multiqc_report_plots
    touch multiqc_report.html
    touch multiqc_report_data/multiqc_stub.txt
    touch multiqc_report_plots/multiqc_stub.txt
    """
}
