process CREATE_TRANSCRIPTOME_FASTA {
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"

    input:
        path(fasta)
        path(gff) 
        tuple val(features), val(attribute)

    output:
        path "${outfile_name}"

    script:
        outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_transcriptome.fasta")
        """
        python $workflow.projectDir/bin/gff_to_fasta_transcriptome.py -fasta $fasta -gff $gff  -f $features -a $attribute -o $outfile_name
        """
}