include { CREATE_TRANSCRIPTOME_FASTA_HOST } from '../../modules/local/create_transcriptome_fasta_host'

include { CREATE_TRANSCRIPTOME_FASTA as CREATE_TRANSCRIPTOME_FASTA_HOST_TRNA } from '../../modules/local/create_transcriptome_fasta'
include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF_TRNA      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF_TRNA_FILE  } from '../../modules/local/uncompress_gff'

include { COMBINE_FILES }  from '../../modules/local/combine_files'


workflow PREPARE_HOST_TRANSCRIPTOME {
  take:
    uncompressed_host_fasta_genome
    ch_gff_host
    ch_gff_host_tRNA

  main:
    
    if(params.read_transcriptome_fasta_host_from_file){
        ch_transcriptome_host        = params.transcriptome_host   ? Channel.fromPath( params.transcriptome_host, checkIfExists: true ) : Channel.empty()
    } else {
        UNCOMPRESS_HOST_GFF(ch_gff_host)

        CREATE_TRANSCRIPTOME_FASTA_HOST(
            uncompressed_host_fasta_genome,
            UNCOMPRESS_HOST_GFF.out
        )        
        if(params.gff_host_tRNA){

            UNCOMPRESS_HOST_GFF_TRNA(UNCOMPRESS_HOST_GFF.out)
            UNCOMPRESS_HOST_GFF_TRNA_FILE(ch_gff_host_tRNA)
            
            parameters = Channel.value(
                [
                    params.gene_feature_gff_to_create_transcriptome_pathogen,
                    params.gene_attribute_gff_to_create_transcriptome_pathogen
                ]
            )
            CREATE_TRANSCRIPTOME_FASTA_HOST_TRNA(
                uncompressed_host_fasta_genome,
                UNCOMPRESS_HOST_GFF_TRNA_FILE.out,
                parameters
            )
            transciptiome_fasta_to_combine = CREATE_TRANSCRIPTOME_FASTA_HOST.out.mix(
                CREATE_TRANSCRIPTOME_FASTA_HOST_TRNA.out
            ).collect()
            COMBINE_FILES(transciptiome_fasta_to_combine)
            ch_transcriptome_host = COMBINE_FILES.out
        } else {
            ch_transcriptome_host = CREATE_TRANSCRIPTOME_FASTA_HOST.out
        }
    }
  emit:
    ch_transcriptome_host
}
