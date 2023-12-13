process SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED {
    tag "split_quants_uniq_mapped_host"
    label 'process_high'
    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
    
    input:
    path quant_table
    path host_annotations
    path pathogen_annotations

    output:
    path 'host_quantification_uniquely_mapped_htseq.tsv', emit: host
    path 'pathogen_quantification_uniquely_mapped_htseq.tsv', emit: pathogen
    env host_tab, emit: scatterplots_host_htseq
    env pathonen_tab, emit: scatterplots_pathogen_htseq
    path "versions.yml"           , emit: versions

    script:
    """
    $workflow.projectDir/bin/split_quant_tables.sh $quant_table $host_annotations $pathogen_annotations quantification_uniquely_mapped_htseq.tsv
        pathonen_tab=\$(if [ \$(cat pathogen_quantification_uniquely_mapped_htseq.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
        host_tab=\$(if [ \$(cat host_quantification_uniquely_mapped_htseq.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}