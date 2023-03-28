include { STAR_GENOMEGENERATE    } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN             } from '../../modules/nf-core/star/align/main'  

workflow SALMON_ALIGNMENT_BASE {

    take:
        ch_reads            // channel: [ val(meta), [ reads ] ]
        ch_genome_fasta     // channel: /path/to/genome.fasta 
        ch_transcript_fasta // channel: /path/to/transcript.fasta
        ch_gtf              // channel: /path/to/genome.gtf 
    main:

        ch_versions = Channel.empty()

        STAR_GENOMEGENERATE ( ch_genome_fasta, ch_gtf )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

        STAR_ALIGN ( ch_reads, STAR_GENOMEGENERATE.out.index, ch_gtf, true, '', '' )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    emit:
        versions = ch_versions                     // channel: [ versions.yml ]
}

