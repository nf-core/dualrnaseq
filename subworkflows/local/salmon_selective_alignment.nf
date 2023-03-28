include { SALMON_INDEX                          } from '../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT                          } from '../../modules/nf-core/salmon/quant/main'
include { COMBINE_QUANTIFICATION_RESULTS_SALMON } from '../../modules/local/combine_quantification_results_salmon/main'
include { SALMON_SPLIT_TABLE                    } from '../../modules/local/salmon_split_table'

workflow SALMON_SELECTIVE_ALIGNMENT {

    take:
        ch_reads            // channel: [ val(meta), [ reads ] ]
        ch_genome_fasta     // channel: /path/to/genome.fasta
        ch_transcript_fasta // channel: /path/to/transcript.fasta
        ch_gtf              // channel: /path/to/genome.gtf
        ch_transcript_fasta_pathogen
        ch_transcript_fasta_host

    main:
        ch_versions = Channel.empty()
        ch_salmon_index = Channel.empty()
        ch_salmon_index = SALMON_INDEX ( ch_genome_fasta, ch_transcript_fasta ).index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

        ch_salmon_quant = Channel.empty()
        def alignment_mode = false
        SALMON_QUANT(ch_reads, ch_salmon_index.collect(), ch_gtf.collect(), ch_transcript_fasta.collect(), alignment_mode, params.libtype)
        
        input_files = SALMON_QUANT.out.results.map{it -> it[1]}.collect()
  //      input_files.view()
        COMBINE_QUANTIFICATION_RESULTS_SALMON(input_files, Channel.value("both"))

        SALMON_SPLIT_TABLE(SALMON_QUANT.out.quant, ch_transcript_fasta_pathogen, ch_transcript_fasta_host)

    emit:
        versions = ch_versions                     // channel: [ versions.yml ]
}
