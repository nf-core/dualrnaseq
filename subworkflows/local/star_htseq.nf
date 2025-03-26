include { HTSEQ } from '../../modules/nf-core/htseq/main'

workflow STAR_HTSEQ {

    take:
        ch_gtf  // channel: /path/to/genome.gtf
    main:

        ch_versions = Channel.empty()
        ch_htseq = Channel.empty()

        HTSEQ ( ch_htseq, ch_gtf )
        ch_versions = ch_versions.mix(HTSEQ.out.versions.first())

    emit:
        versions = ch_versions  // channel: [ versions.yml ]
}

