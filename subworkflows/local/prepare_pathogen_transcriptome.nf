include { CREATE_TRANSCRIPTOME_FASTA as CREATE_TRANSCRIPTOME_FASTA_PATHOGEN } from '../../modules/local/create_transcriptome_fasta'
include { COMBINE_FILES } from '../../modules/local/combine_files'


workflow PREPARE_PATHOGEN_TRANSCRIPTOME {
    take:
        uncompress_pathogen_fasta_genome
        uncompress_pathogen_gff
    main:
        if(params.read_transcriptome_fasta_pathogen_from_file){
            transcriptome_pathogen  = params.transcriptome_pathogen ? Channel.fromPath( params.transcriptome_pathogen, checkIfExists: true ) : Channel.empty()
        } else {
            parameters = Channel.value(
                [
                    params.gene_feature_gff_to_create_transcriptome_pathogen,
                    params.gene_attribute_gff_to_create_transcriptome_pathogen
                ]
            )
            CREATE_TRANSCRIPTOME_FASTA_PATHOGEN(
                uncompress_pathogen_fasta_genome,
                uncompress_pathogen_gff,
                parameters
            )
            transcriptome_pathogen = CREATE_TRANSCRIPTOME_FASTA_PATHOGEN.out
        }
    emit:
        transcriptome_pathogen

}