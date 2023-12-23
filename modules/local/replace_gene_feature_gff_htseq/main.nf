process REPLACE_GENE_FEATURE_GFF_HOST_HTSEQ {
    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
    storeDir "${params.outdir}/references"
    tag "repl_gene_feature_gff_host"

    label 'process_high'

    input:
    path(gff)
    val(features)

    output:
    path "${outfile_name}" into combine_gff_host 

    script:
    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_quant_feature.gff3")
    """
    $workflow.projectDir/bin/replace_feature_gff.sh $gff ${outfile_name} $features

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
}
