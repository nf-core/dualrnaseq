process PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH{
    publishDir "${params.outdir}/mapping_statistics/${tool}/RNA_classes_pathogen", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/${tool}/RNA_classes_pathogen"
    tag "plot_rna_stats_htseq_pathogen"

    label 'process_high'

    input:
    file stats_table from plot_RNA_stats_pathogen_htseq_u_m
    val plot_rna from plot_RNA_stats_pathogen_htseq_u_m_boolean
    val(tool)

    output:
    file "*.pdf"

    when:
    plot_rna.toBoolean()

    script:
    """
    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py -i $stats_table
    """
}

		