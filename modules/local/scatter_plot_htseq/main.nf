process SCATTER_PLOT_HTSEQ {
    publishDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots"
    tag "scatter_plot_pathogen_htseq"

    label 'process_high'

    input:
    file quant_table from quant_scatter_plot_pathogen_htseq_u_m
    val attribute from atr_scatter_plot_pathogen_htseq_u_m
    val replicates from repl_scatter_plots_htseq_pathogen
    val pathogen_table_non_empty from scatterplots_pathogen_htseq

    output:
    file ('*.pdf')

    when:
    replicates.toBoolean()
    pathogen_table_non_empty.toBoolean()

    script:
    """
    python $workflow.projectDir/bin/scatter_plots.py -q $quant_table -a $attribute -org pathogen
    """
}