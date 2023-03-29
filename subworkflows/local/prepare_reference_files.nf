include { UNCOMPRESS_FASTA_GENOME as UNCOMPRESS_HOST_FASTA_GENOME      } from '../../modules/local/uncompress_fasta_genome'
include { UNCOMPRESS_FASTA_GENOME as UNCOMPRESS_PATHOGEN_FASTA_GENOME  } from '../../modules/local/uncompress_fasta_genome'

include { PREPARE_HOST_TRANSCRIPTOME      } from '../../subworkflows/local/prepare_host_transcriptome'
include { PREPARE_PATHOGEN_TRANSCRIPTOME  } from '../../subworkflows/local/prepare_pathogen_transcriptome'
include { COMBINE_FILES as COMBINE_FILES_GENOMES } from '../../modules/local/combine_files'

include { COMBINE_FILES as COMBINE_FILES_TRANSCRIPTOME_FILES } from '../../modules/local/combine_files'


workflow PREPARE_REFERENCE_FILES{
  take: 
    ch_fasta_host
    ch_gff_host
    ch_gff_host_tRNA
    ch_fasta_pathogen
    ch_gff_pathogen

  main:
    ch_transcriptome = Channel.empty()
    ch_host_transcriptome = Channel.empty()
    ch_pathogen_transcriptome = Channel.empty()

    // uncompress genome
    UNCOMPRESS_PATHOGEN_FASTA_GENOME(ch_fasta_pathogen)
    UNCOMPRESS_HOST_FASTA_GENOME(ch_fasta_host)
    ch_genome_fasta = UNCOMPRESS_HOST_FASTA_GENOME.out

    // combine pathogen and host genome
    fasta_files_to_combine = UNCOMPRESS_HOST_FASTA_GENOME.out.mix(
        UNCOMPRESS_PATHOGEN_FASTA_GENOME.out
    ).collect()
    COMBINE_FILES_GENOMES(fasta_files_to_combine)    

    if(params.run_salmon_selective_alignment | params.run_salmon_alignment_based_mode) {

      PREPARE_HOST_TRANSCRIPTOME(
        UNCOMPRESS_HOST_FASTA_GENOME.out,
        ch_gff_host,
        ch_gff_host_tRNA
      )
      ch_host_transcriptome = PREPARE_HOST_TRANSCRIPTOME.out

      PREPARE_PATHOGEN_TRANSCRIPTOME(
        UNCOMPRESS_PATHOGEN_FASTA_GENOME.out.collect(),
        ch_gff_pathogen
      )
      ch_pathogen_transcriptome = PREPARE_PATHOGEN_TRANSCRIPTOME.out

      // combine pathogen and host transcriptome
      transciptiome_transcriptome_to_combine = ch_host_transcriptome.mix(
        ch_pathogen_transcriptome
      ).collect()
      COMBINE_FILES_TRANSCRIPTOME_FILES(
        transciptiome_transcriptome_to_combine
      )
      ch_transcriptome = COMBINE_FILES_TRANSCRIPTOME_FILES.out
    }
    emit:
      genome_fasta = ch_genome_fasta
      transcript_fasta = ch_transcriptome
      transcript_fasta_host = ch_host_transcriptome
      transcript_fasta_pathogen = ch_pathogen_transcriptome
}
