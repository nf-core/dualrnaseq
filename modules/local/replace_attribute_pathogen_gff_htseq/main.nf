process REPLACE_ATTRIBUTE_PATHOGEN_GFF_HTSEQ {
    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
    storeDir "${params.outdir}/references"
    tag "repl_attribute_pathogen_gff"

    label 'process_high'

    input:
    path(gff)
    val(host_attribute)
    val(pathogen_attribute)

    output:
    path "${outfile_name}", emit: gff3

    script:
    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_new_attribute.gff3")
    """
    $workflow.projectDir/bin/replace_attribute_gff.sh $gff ${outfile_name} $host_attribute $pathogen_attribute

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
}