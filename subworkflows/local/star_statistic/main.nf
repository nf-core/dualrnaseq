include { REMOVE_CROSSMAPPED_READS             } from '../../../modules/local/remove_crossmapped_reads/main'
include { EXTRACT_PROCESSED_READS              } from '../../../modules/local/extract_processed_reads'
include { COLLECT_PROCESSED_READS_STAR         } from '../../../modules/local/collect_processed_reads_star/main'
include { UNIQUE_MAPPING_STATS_STAR            } from '../../../modules/local/unique_mapping_stats_star/main'
include { COLLECT_STATS_STAR_UNIQUELY_MAPPED   } from '../../../modules/local/collect_stats_star_uniquely_mapped/main'
include { COUNT_CROSSMAPPED_READS              } from '../../../modules/local/count_crossmapped_reads/main'

include { MULTI_MAPPING_STATS                  } from '../../../modules/local/multi_mapping_stats/main'
include { STAR_MAPPING_STATS                  } from '../../../modules/local/star_mapping_stats/main'
include { COLLECT_STATS_STAR_MULTI_MAPPED      } from '../../../modules/local/collect_stats_star_multi_mapped/main'
include { PLOT_STAR_MAPPING_STATS }  from '../../../modules/local/plot_star_mapping_stats/main' 
workflow STAR_STATISTIC {
    take:
        star_bam
        star_log_final
        reference_host_name
        reference_pathogen_name
        total_read_pairs
    main:
        ch_versions = Channel.empty()
        
        if(params.mapping_statistics) {

            REMOVE_CROSSMAPPED_READS(
                star_bam, 
                reference_host_name,
                reference_pathogen_name
            )
            ch_versions = ch_versions.mix(REMOVE_CROSSMAPPED_READS.out.versions.first())
            
            EXTRACT_PROCESSED_READS(star_log_final, "star")
            ch_versions = ch_versions.mix(EXTRACT_PROCESSED_READS.out.versions.first())

            COLLECT_PROCESSED_READS_STAR(EXTRACT_PROCESSED_READS.out.collect_results.collect())
            ch_versions = ch_versions.mix(COLLECT_PROCESSED_READS_STAR.out.versions.first())
            
            UNIQUE_MAPPING_STATS_STAR(
                star_bam,         
                reference_host_name, 
                reference_pathogen_name
            )
            ch_versions = ch_versions.mix(UNIQUE_MAPPING_STATS_STAR.out.versions.first())

            COLLECT_STATS_STAR_UNIQUELY_MAPPED(UNIQUE_MAPPING_STATS_STAR.out.txt)
            ch_versions = ch_versions.mix(COLLECT_STATS_STAR_UNIQUELY_MAPPED.out.versions.first())

            COUNT_CROSSMAPPED_READS(REMOVE_CROSSMAPPED_READS.out.count_crossmapped_reads)
            ch_versions = ch_versions.mix(COUNT_CROSSMAPPED_READS.out.versions.first())

            MULTI_MAPPING_STATS(
                REMOVE_CROSSMAPPED_READS.out.alignment_multi_mapping_stats, 
                reference_host_name, 
                reference_pathogen_name
            )
            ch_versions = ch_versions.mix(MULTI_MAPPING_STATS.out.versions.first())

            COLLECT_STATS_STAR_MULTI_MAPPED(MULTI_MAPPING_STATS.out.txt)
            ch_versions = ch_versions.mix(COLLECT_STATS_STAR_MULTI_MAPPED.out.versions.first())

            STAR_MAPPING_STATS(
                total_read_pairs,
                COLLECT_PROCESSED_READS_STAR.out.tsv,
                COLLECT_STATS_STAR_UNIQUELY_MAPPED.out.tsv,
                COLLECT_STATS_STAR_MULTI_MAPPED.out.tsv,
                COUNT_CROSSMAPPED_READS.out.txt
            )
            ch_versions = ch_versions.mix(STAR_MAPPING_STATS.out.versions.first())

            PLOT_STAR_MAPPING_STATS(STAR_MAPPING_STATS.out.tsv)
            ch_versions = ch_versions.mix(PLOT_STAR_MAPPING_STATS.out.versions.first())
        }
    emit:
        versions = ch_versions  // channel: [ versions.yml ]
}
