include { HTSEQ_COUNT } from '../../../modules/local/htseq'
include { 
    COMBINE_ANNOTATIONS_QUANT_UNIQUELY_MAPPED_HOST as COMBINE_ANNOTATIONS_QUANT_PATHOGEN_UNIQUELY_MAPPED_HOST;
    COMBINE_ANNOTATIONS_QUANT_UNIQUELY_MAPPED_HOST as COMBINE_ANNOTATIONS_QUANT_HOST_UNIQUELY_MAPPED_HOST; 
} from  '../../../modules/local/combine_annotations_quant_uniquely_mapped_host/main'

include { COMBINE_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED } from  '../../../modules/local/combine_quantification_tables_htseq_uniquely_mapped/main'


include { HTSEQ_QUANTIFICATION_STATS_UNIQUELY_MAPPED } from  '../../../modules/local/htseq_quantification_stats_uniquely_mapped/main'
include { HTSEQ_UNIQUELY_MAPPED_TPM } from  '../../../modules/local/htseq_uniquely_mapped_tpm/main'

include { SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED } from  '../../../modules/local/split_quantification_tables_htseq_uniquely_mapped/main'
 
include { PLOT_MAPPING_STATS_HOST_PATHOGEN_HTSEQ_UNIQUELY_MAPPED } from  '../../../modules/local/plot_mapping_stats_host_pathogen_htseq_uniquely_mapped/main'
include { RNA_STATISTICS } from '../rna_statistics/main'


workflow STAR_HTSEQ {

    take:
        star_bam // STAR_ALIGN.out.bam_unsorted, 
        quantification // PREPARE_REFERENCE_FILES.quantification_gff_u_m,
        // annotations_host_htseq
        // annotations_pathogen_htseq
        gff_host // PREPARE_REFERENCE_FILES.gff_host,
        gff_pathogen // PREPARE_REFERENCE_FILES.gff_pathogen,
        // star_mapping_stats_tsv

    main:
        ch_versions = Channel.empty()
        ch_htseq = Channel.empty()

        HTSEQ_COUNT ( star_bam, quantification, "quant" )
        ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())

	    COMBINE_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED(
            HTSEQ_COUNT.out.txt.map {meta, path -> path }.collect(), params.host_gff_attribute
        )
        ch_versions = ch_versions.mix(COMBINE_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.versions.first())

        HTSEQ_UNIQUELY_MAPPED_TPM(
            COMBINE_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.tsv,
            params.host_gff_attribute,
            gff_host,
            gff_pathogen
        )
        ch_versions = ch_versions.mix(HTSEQ_UNIQUELY_MAPPED_TPM.out.versions.first())

        // SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED(
        //     HTSEQ_UNIQUELY_MAPPED_TPM.out.tsv,
        //     annotations_host_htseq,
        //     annotations_pathogen_htseq
        // )
        // ch_versions = ch_versions.mix(SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.versions.first())

        // COMBINE_ANNOTATIONS_QUANT_PATHOGEN_UNIQUELY_MAPPED_HOST(
        //     SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.pathogen,
        //     annotations_pathogen_htseq,
        //     params.host_gff_attribute,
        //     "pathogen"
        // )
        
        // COMBINE_ANNOTATIONS_QUANT_HOST_UNIQUELY_MAPPED_HOST(
        //    SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.host,
        //     annotations_host_htseq,
        //     params.host_gff_attribute,
        //     "host"
        // )
        // if(params.mapping_statistics) {
            
            

        //     HTSEQ_QUANTIFICATION_STATS_UNIQUELY_MAPPED(
        //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.host,
        //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.pathogen,
        //         Channel.value(params.host_gff_attribute),
        //         star_mapping_stats_tsv
        //     )
        //     ch_versions = ch_versions.mix(HTSEQ_QUANTIFICATION_STATS_UNIQUELY_MAPPED.out.versions.first())

        //     PLOT_MAPPING_STATS_HOST_PATHOGEN_HTSEQ_UNIQUELY_MAPPED(
        //         HTSEQ_QUANTIFICATION_STATS_UNIQUELY_MAPPED.out.tsv
        //     )
        //     ch_versions = ch_versions.mix(PLOT_MAPPING_STATS_HOST_PATHOGEN_HTSEQ_UNIQUELY_MAPPED.out.versions.first())
        //     RNA_STATISTICS(
        //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.host,
        //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.pathogen,
        //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.scatterplots_host_htseq,
        //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.scatterplots_pathogen_htseq,
        //         annotations_host_htseq,
        //         annotations_pathogen_htseq,
        //         Channel.value(params.host_gff_attribute),
        //         "htseq"
        //     )
        //     ch_versions = ch_versions.mix(RNA_STATISTICS.out.versions.first())

        // }
    emit:
        versions = ch_versions  // channel: [ versions.yml ]
}

workflow {
    STAR_HTSEQ()
}