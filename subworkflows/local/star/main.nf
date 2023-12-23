include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN                        } from '../../../modules/nf-core/star/align/main'

workflow STAR {
    take:
        ch_reads
        ch_genome_fasta
        ch_host_gff
    main:
        ch_versions = Channel.empty()
        STAR_GENOMEGENERATE ( ch_genome_fasta, ch_host_gff )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

        STAR_ALIGN ( ch_reads, STAR_GENOMEGENERATE.out.index, ch_host_gff, true, '', '' )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
    emit:
        bam = STAR_ALIGN.out.bam
        log_final = STAR_ALIGN.out.log_final
        versions = ch_versions
}