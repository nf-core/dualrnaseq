include { CREATE_TRANSCRIPTOME_FASTA } from '../../modules/local/create_transcriptome_fasta'
include { COMBINE_FILES } from '../../modules/local/combine_files'
include { UNCOMPRESS_GFF} from '../../modules/local/uncompress_gff'

include { REPLACE_GENE_FEATURE_GFF_SALMON } from '../../modules/local/replace_gene_feature'

workflow PREPARE_PATHOGEN_TRANSCRIPTOME {
    take:
        uncompress_fasta_genome
        ch_gff_pathogen
    main:
        UNCOMPRESS_GFF (ch_gff_pathogen)
        ch_uncompress_gff = UNCOMPRESS_GFF.out

        if(params.read_transcriptome_fasta_pathogen_from_file){
            ch_transcriptome  = params.transcriptome_pathogen ? Channel.fromPath( params.transcriptome_pathogen, checkIfExists: true ) : Channel.empty()
        } else {
            parameters = Channel.value(
                [
                    params.gene_feature_gff_to_create_transcriptome_pathogen,
                    params.gene_attribute_gff_to_create_transcriptome_pathogen
                ]
            )
            CREATE_TRANSCRIPTOME_FASTA(
                uncompress_fasta_genome,
                ch_uncompress_gff,
                parameters
            )

            ch_transcriptome = CREATE_TRANSCRIPTOME_FASTA.out
        }

    emit:
        transcriptome = ch_transcriptome
        uncompress_gff = ch_uncompress_gff

}