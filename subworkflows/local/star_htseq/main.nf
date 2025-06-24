include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate'
include { STAR_ALIGN          } from '../../../modules/local/star_align_genome'
include { HTSEQ_COUNT         } from '../../../modules/local/htseq_count'

workflow STAR_HTSEQ {
    take:
    ch_reads                      // channel: [ val(meta), [ reads ] ]
    ch_host_pathogen_fasta_genome
    ch_host_pathogen_gff

    main:

    ch_versions = Channel.empty()

    // -------
    // Run create STAR index
    // -------
    STAR_GENOMEGENERATE(
        ch_host_pathogen_fasta_genome,
        ch_host_pathogen_gff,
    )
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)


    // -------
    // Run STAR align
    // -------
    STAR_ALIGN(
        ch_reads,
        STAR_GENOMEGENERATE.out.index,
        ch_host_pathogen_gff,
        true,
        '',
        '',
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)


    // -------
    // Run HTSeq-count
    // -------
    if (params.htseq) {
        HTSEQ_COUNT(
            STAR_ALIGN.out.bam_sorted,
            ch_host_pathogen_gff,
        )
        ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())
    }

    // -------
    // Split each count table into host and pathogen reads (for each dataset)
    // -------
    // SALMON_SPLIT_TABLE_EACH(SALMON_QUANT.out.quant, //HTSEQ_COUNT.out.counts
    //                         ch_pathogen_fasta_transcripts,
    //                         ch_host_fasta_transcripts
    //                         )

    // -------
    //  Combine count results from all samples
    // -------

    // -------
    //  Separate out host and pathogen reads from combined counts
    // -------

    // -------
    // Combine all meta data from each datasets
    // -------



    // -------
    //  Capture the number of counted reads by HTSeq and save as output
    // -------

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
