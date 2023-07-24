process CREATE_TRANSCRIPTOME_FASTA_GFFREAD {
    tag "create_transcripts_host"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/nfcore/dualrnaseq:dev' :
        'docker.io/nfcore/dualrnaseq:dev' }"

    input:
        path(fasta)
        path(gff)
    output:
        path "${outfile_name}"

    script:
    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_transcriptome.fasta")
    """
    gffread -w $outfile_name -g $fasta $gff
    """
}
