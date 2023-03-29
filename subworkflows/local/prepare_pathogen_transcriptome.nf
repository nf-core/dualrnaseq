include { CREATE_TRANSCRIPTOME_FASTA as CREATE_TRANSCRIPTOME_FASTA_PATHOGEN } from '../../modules/local/create_transcriptome_fasta'
include { COMBINE_FILES } from '../../modules/local/combine_files'
include { UNCOMPRESS_GFF as UNCOMPRESS_PATHOGEN_GFF  } from '../../modules/local/uncompress_gff'

include { EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_PATHOGEN_HTSEQ } from '../../modules/local/extract_annotations'

workflow PREPARE_PATHOGEN_TRANSCRIPTOME {
    take:
        uncompress_pathogen_fasta_genome
        ch_gff_pathogen
    main:
        if(params.read_transcriptome_fasta_pathogen_from_file){
            transcriptome_pathogen  = params.transcriptome_pathogen ? Channel.fromPath( params.transcriptome_pathogen, checkIfExists: true ) : Channel.empty()
        } else {
            UNCOMPRESS_PATHOGEN_GFF(ch_gff_pathogen)
            
            EXTRACT_ANNOTATIONS_PATHOGEN_HTSEQ (
                UNCOMPRESS_PATHOGEN_GFF.out,
                params.gene_feature_gff_to_quantify_pathogen,
                params.pathogen_gff_attribute,
                params.extract_annotations_pathogen_htseq_organism,
                'htseq'
            )
            parameters = Channel.value(
                [
                    params.gene_feature_gff_to_create_transcriptome_pathogen,
                    params.gene_attribute_gff_to_create_transcriptome_pathogen
                ]
            )
            CREATE_TRANSCRIPTOME_FASTA_PATHOGEN(
                uncompress_pathogen_fasta_genome,
                UNCOMPRESS_PATHOGEN_GFF.out,
                parameters
            )
            transcriptome_pathogen = CREATE_TRANSCRIPTOME_FASTA_PATHOGEN.out
        }



    emit:
        transcriptome_pathogen

}