process EXTRA_REFERENCE_NAME_STAR {
    tag "extract_ref_names_host_star"
    label 'process_high'

    publishDir "${params.outdir}/references", mode: params.publish_dir_mode 
    storeDir "${params.outdir}/references"

    input:
    path(fasta) from genome_fasta

    output:
    file "reference_names.txt" into reference_name

    script:
    """
    $workflow.projectDir/bin/extract_reference_names_from_fasta_files.sh reference_names.txt $fasta
    """
}