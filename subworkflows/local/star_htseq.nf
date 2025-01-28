include { STAR_GENOMEGENERATE                               } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN                                        } from '../../modules/local/star_align_genome'
include { HTSEQ_COUNT                                       } from '../../modules/local/htseq_count'

workflow STAR_HTSEQ {

    take:
        ch_reads            // channel: [ val(meta), [ reads ] ]
        ch_host_pathogen_fasta_genome
        ch_host_pathogen_gff
    main:

        ch_versions = Channel.empty()

        // -------
        // Run create STAR index
        // -------
        STAR_GENOMEGENERATE (
            ch_host_pathogen_fasta_genome,
            ch_host_pathogen_gff
            )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)


        // -------
        // Run STAR align
        // -------
        STAR_ALIGN ( ch_reads, // reads
                     STAR_GENOMEGENERATE.out.index, // index
                     ch_host_pathogen_gff, // GTF
                     true, //star_ignore_sjdbgtf
                     '', // seq_platform
                     '' // seq_centre
                    )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)


        // -------
        // Run HTSeq-count
        // -------
        // if ( params.run_htseq ) {

        //     HTSEQ_COUNT (
        //         STAR_ALIGN.out.bam_sorted

        //         ch_gtf
        //         )
        //     ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())
        // }



    emit:
    versions = ch_versions  // channel: [ versions.yml ]
}

