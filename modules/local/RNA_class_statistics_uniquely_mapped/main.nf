process RNA_CLASS_STATISTICS_UNIQUELY_MAPPED {
    publishDir "${params.outdir}/mapping_statistics/${tool}/RNA_classes_host", mode: params.publish_dir_mode
    tag "rna_class_stats_htseq_host"

    label 'process_high'

    input:
    path(quant_table)
    val(attribute)
    path(gene_annotations)
    path(rna_classes_to_replace)
    val(organism)
    val(tool)

    output:
    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host_htseq_u_m
    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host_combined_htseq_u_m
    file "host_RNA_classes_sum_counts_*.tsv"
    file "host_gene_types_groups_*"
    stdout plot_RNA_stats_host_htseq_u_m_boolean
    stdout plot_RNA_stats_host_combined_htseq_u_m_boolean

    shell:
    '''
    python !{workflow.projectDir}/bin/RNA_class_content.py \\
    -q !{quant_table} \\
    -a !{attribute} \\
    -annotations !{gene_annotations} \\
    -rna !{rna_classes_to_replace} \\
    -q_tool !{tool} \\
    -org !{organism} 2>&1
    '''
}
