include { SALMON_INDEX      } from '../../modules/nf-core/salmon/index/main'
include { SALMON_QUANT     } from '../../modules/nf-core/salmon/quant/main'

workflow SALMON_SELECTIVE_ALIGNMENT {

    take:
    ch_reads            // channel: [ val(meta), [ reads ] ]
    ch_genome_fasta     // channel: /path/to/genome.fasta
    ch_transcript_fasta // channel: /path/to/transcript.fasta
    ch_gtf              // channel: /path/to/genome.gtf

    main:
    ch_versions = Channel.empty()
    ch_salmon_index = Channel.empty()
    ch_salmon_index = SALMON_INDEX ( ch_genome_fasta, ch_transcript_fasta ).index
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    ch_salmon_quant = Channel.empty()
    def alignment_mode = false
    SALMON_QUANT(ch_reads, ch_salmon_index, ch_gtf, ch_transcript_fasta, alignment_mode, params.libtype)

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
