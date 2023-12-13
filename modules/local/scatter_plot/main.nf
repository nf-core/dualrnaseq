process SCATTER_PLOT {
    publishDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots"
    tag "scatter_plot_pathogen_htseq"

    label 'process_high'

    input:
    path quant_table   
    val attribute
    val replicates
    val pathogen_table_non_empty
    val host

    output:
    path ('*.pdf')
    path "versions.yml"           , emit: versions

    when:
    replicates.toBoolean()
    pathogen_table_non_empty.toBoolean()

    script:
    """
    python $workflow.projectDir/bin/scatter_plots.py \\
    -q $quant_table \\
    -a $attribute \\
    -org $host
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}