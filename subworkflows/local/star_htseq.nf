include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN                        } from '../../modules/nf-core/star/align/main'  
include { HTSEQ                             } from '../../modules/local/htseq'
include { COMBINE_QUANTIFICATION_RESULTS_HTS} from '../../modules/local/combine_quantification_results_hts'
include { HTSEQ_UNIQUELY_MAPPED_TPM         } from '../../modules/local/htseq_uniquely_mapped_tpm'
include { SPLIT_QUATIFICATION_TABLES_HTSEQ  } from '../../modules/local/split_quantification_tables_htseq_uniquely_mapped'
include { 
            COMBINE_ANNOTATIONS_QUANT as COMBINE_ANNOTATIONS_QUANT_PATHOGEN;
            COMBINE_ANNOTATIONS_QUANT as COMBINE_ANNOTATIONS_QUANT_HOST
        } from '../../modules/local/combine_annotations_quant'


workflow STAR_HTSEQ {

    take:
        ch_reads            // channel: [ val(meta), [ reads ] ]
        ch_genome_fasta     // channel: /path/to/genome.fasta 
        ch_transcript_fasta // channel: /path/to/transcript.fasta
        ch_gff_host_pathogen              // channel: /path/to/genome.gtf 
        ch_transcript_fasta_pathogen
        ch_transcript_fasta_host
        ch_gff_host
        ch_gff_pathogen
        ch_gff_host_unzipped
        annotations_host_htseq
        annotations_pathogen_htseq

    main:
        ch_versions = Channel.empty()

        STAR_GENOMEGENERATE ( ch_genome_fasta, ch_gff_host_unzipped )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

        STAR_ALIGN ( ch_reads, STAR_GENOMEGENERATE.out.index, ch_gff_host_unzipped, true, '', '' )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

        HTSEQ ( STAR_ALIGN.out.bam, 
                ch_gff_host_pathogen, 
                params.host_gff_attribute,
                params.htseq_stranded
                )
        ch_versions = ch_versions.mix(HTSEQ.out.versions.first())

        input_files = HTSEQ.out.results.map{it -> it[1]}.collect()
        COMBINE_QUANTIFICATION_RESULTS_HTS(input_files, params.host_gff_attribute)

        COMBINE_QUANTIFICATION_RESULTS_HTS.out.combined_quant_data
        .map {it ->
            def meta = [:]
            meta.id  = "combined"
            path_res = it
            return [ meta, [ it ] ]
        }.set{ combined_htseq_quant }


        HTSEQ_UNIQUELY_MAPPED_TPM(
            combined_htseq_quant,
            params.host_gff_attribute,
            ch_gff_host,
            ch_gff_pathogen
        )

        SPLIT_QUATIFICATION_TABLES_HTSEQ(
            HTSEQ_UNIQUELY_MAPPED_TPM.out,
            annotations_host_htseq,
            annotations_pathogen_htseq
        )

        COMBINE_ANNOTATIONS_QUANT_PATHOGEN(
            SPLIT_QUATIFICATION_TABLES_HTSEQ.out.pathogen_quantification_stats_htseq,
            annotations_pathogen_htseq,
            params.host_gff_attribute
        )

        COMBINE_ANNOTATIONS_QUANT_HOST(
            SPLIT_QUATIFICATION_TABLES_HTSEQ.out.host_quantification_stats_htseq,
            annotations_host_htseq,
            params.host_gff_attribute
        )

    emit:
        versions = ch_versions  // channel: [ versions.yml ]
}

