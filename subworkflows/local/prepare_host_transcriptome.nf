include { CREATE_TRANSCRIPTOME_FASTA_HOST } from '../../modules/local/create_transcriptome_fasta_host'

include { CREATE_TRANSCRIPTOME_FASTA as CREATE_TRANSCRIPTOME_FASTA_HOST_TRNA } from '../../modules/local/create_transcriptome_fasta'

include { COMBINE_FILES }  from '../../modules/local/combine_files'


workflow PREPARE_HOST_TRANSCRIPTOME {
  take:
    uncompressed_host_fasta_genome
    uncompressed_host_gff

  main:
    if(params.read_transcriptome_fasta_host_from_file){
        ch_transcriptome_host        = params.transcriptome_host   ? Channel.fromPath( params.transcriptome_host, checkIfExists: true ) : Channel.empty()
    } else {
        CREATE_TRANSCRIPTOME_FASTA_HOST(
            uncompressed_host_fasta_genome,
            uncompressed_host_gff
            
        )        
        if(params.gff_host_tRNA){

            UNCOMPRESS_HOST_GFF_TRNA(ch_gff_host)
            UNCOMPRESS_HOST_GFF_TRNA_FILE(ch_gff_host_tRNA)
            
            parameters = Channel.value(
                [
                    params.gene_feature_gff_to_create_transcriptome_pathogen,
                    params.gene_attribute_gff_to_create_transcriptome_pathogen
                ]
            )
            CREATE_TRANSCRIPTOME_FASTA_HOST_TRNA(
                UNCOMPRESS_HOST_FASTA_GENOME.out,
                UNCOMPRESS_HOST_GFF_TRNA_FILE.out,
                parameters
            )
            transciptiome_fasta_to_combine = CREATE_TRANSCRIPTIOME_FASTA_HOST.out.mix(
                CREATE_TRANSCRIPTIOME_FASTA_HOST_TRNA.out
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
