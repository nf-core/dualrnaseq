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
    path "host_RNA_classes_percentage_*.tsv", emit: tsv
    path "host_gene_types_groups_*", emit: random
    path "versions.yml" , emit: versions

    shell:
    '''
    python !{workflow.projectDir}/bin/RNA_class_content.py \\
    -q !{quant_table} \\
    -a !{attribute} \\
    -annotations !{gene_annotations} \\
    -rna !{rna_classes_to_replace} \\
    -q_tool !{tool} \\
    -org !{organism} 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    '''
}
