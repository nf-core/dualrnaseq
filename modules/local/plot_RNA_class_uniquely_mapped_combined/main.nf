process PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED {
    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host"
    tag "plt_rna_stats_htseq_host_all"

    label 'process_high'

    input:
    file stats_table from plot_RNA_stats_host_combined_htseq_u_m
    val plot_rna from plot_RNA_stats_host_combined_htseq_u_m_boolean

    output:
    file "RNA_class_stats_combined_host.pdf"

    when:
    plot_rna.toBoolean()

    script:
    """
    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py -i $stats_table -org host
    """
}
