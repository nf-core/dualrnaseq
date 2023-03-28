include { UNCOMPRESS_FASTA_GENOME as UNCOMPRESS_HOST_FASTA_GENOME      } from '../../modules/local/uncompress_fasta_genome'
include { UNCOMPRESS_FASTA_GENOME as UNCOMPRESS_PATHOGEN_FASTA_GENOME  } from '../../modules/local/uncompress_fasta_genome'

include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_PATHOGEN_GFF  } from '../../modules/local/uncompress_gff'

include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF_TRNA      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF_TRNA_FILE  } from '../../modules/local/uncompress_gff'

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
      UNCOMPRESS_PATHOGEN_FASTA_GENOME(ch_fasta_pathogen)
      UNCOMPRESS_HOST_FASTA_GENOME(ch_fasta_host)


      // combine pathogen and host genome
      fasta_files_to_combine = UNCOMPRESS_HOST_FASTA_GENOME.out.mix(
          UNCOMPRESS_PATHOGEN_FASTA_GENOME.out
      ).collect()
      COMBINE_FILES_GENOMES(fasta_files_to_combine)

      UNCOMPRESS_HOST_GFF(ch_gff_host)
      UNCOMPRESS_PATHOGEN_GFF(ch_gff_pathogen)

      PREPARE_HOST_TRANSCRIPTOME(
        UNCOMPRESS_HOST_FASTA_GENOME.out,
        UNCOMPRESS_HOST_GFF.out
      )

      PREPARE_PATHOGEN_TRANSCRIPTOME(
        UNCOMPRESS_PATHOGEN_FASTA_GENOME.out.collect(),
        UNCOMPRESS_PATHOGEN_GFF.out
      )

      // combine pathogen and host transcriptome
      transciptiome_transcriptome_to_combine = PREPARE_HOST_TRANSCRIPTOME.out.mix(
        PREPARE_PATHOGEN_TRANSCRIPTOME.out
      ).collect()
      COMBINE_FILES_TRANSCRIPTOME_FILES(
        transciptiome_transcriptome_to_combine
      )
}