include { CREATE_TRANSCRIPTOME_FASTA_GFFREAD } from '../../modules/local/create_transcriptome_fasta_gffread'

include { CREATE_TRANSCRIPTOME_FASTA } from '../../modules/local/create_transcriptome_fasta'

include { UNCOMPRESS_GFF as UNCOMPRESS_GFF      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_GFF_TRNA      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_GFF_TRNA_FILE  } from '../../modules/local/uncompress_gff'

include { COMBINE_FILES as COMBINE_FILES_1; COMBINE_FILES as COMBINE_FILES_2 }  from '../../modules/local/combine_files'

include { REPLACE_GENE_FEATURE_GFF_SALMON } from '../../modules/local/replace_gene_feature'



workflow PREPARE_HOST_TRANSCRIPTOME {
  take:
    uncompressed_fasta_genome
    ch_gff_host
    ch_gff_host_tRNA

  main:
    
    if (params.read_transcriptome_fasta_host_from_file) {
        ch_transcriptome_host        = params.transcriptome_host   ? Channel.fromPath( params.transcriptome_host, checkIfExists: true ) : Channel.empty()
    } else {
        UNCOMPRESS_GFF(ch_gff_host)
        CREATE_TRANSCRIPTOME_FASTA_GFFREAD(
            uncompressed_fasta_genome,
            UNCOMPRESS_GFF.out
        )  

        if(params.gff_host_tRNA){
            UNCOMPRESS_GFF_TRNA(ch_gff_host)
            UNCOMPRESS_GFF_TRNA_FILE(ch_gff_host_tRNA)
            parameters = Channel.value(
                [
                    params.gene_feature_gff_to_create_transcriptome_pathogen,
                    params.gene_attribute_gff_to_create_transcriptome_pathogen
                ]
            )
            CREATE_TRANSCRIPTOME_FASTA(
                uncompressed_fasta_genome,
                UNCOMPRESS_GFF_TRNA_FILE.out,
                parameters
            )
            transciptiome_fasta_to_combine = CREATE_TRANSCRIPTOME_FASTA_GFFREAD.out.mix(
                CREATE_TRANSCRIPTOME_FASTA.out
            ).collect()

            COMBINE_FILES_1(transciptiome_fasta_to_combine)

            COMBINE_FILES_2(
                UNCOMPRESS_GFF_TRNA.out.mix(UNCOMPRESS_GFF_TRNA_FILE.out).collect()
            )

            ch_transcriptome = COMBINE_FILES_1.out
            ch_uncompress_gff_trna_file = UNCOMPRESS_GFF_TRNA_FILE.out
            ch_uncompress_gff_trna = UNCOMPRESS_GFF_TRNA.out
            ch_uncompress_gff = Channel.empty()
        } else {
            ch_uncompress_gff_trna_file = Channel.empty()
            ch_uncompress_gff_trna = Channel.empty()
            ch_uncompress_gff = UNCOMPRESS_GFF.out
            ch_transcriptome = CREATE_TRANSCRIPTOME_FASTA_GFFREAD.out
        }
    }
  emit:
    transcriptome = ch_transcriptome
    uncompress_gff_trna_file = ch_uncompress_gff_trna_file
    uncompress_gff_trna = ch_uncompress_gff_trna
    uncompress_gff = ch_uncompress_gff
}
