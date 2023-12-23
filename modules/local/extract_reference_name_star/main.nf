process EXTRA_REFERENCE_NAME_STAR {
    tag "extract_ref_names_host_star"
    label 'process_high'

    publishDir "${params.outdir}/references", mode: params.publish_dir_mode 

    input:
    path(fasta)
    val(organism)

    output:
    path "reference_names_*.txt", emit: txt

    script:
    """
    $workflow.projectDir/bin/extract_reference_names_from_fasta_files.sh reference_names_${organism}.txt $fasta
    """
}