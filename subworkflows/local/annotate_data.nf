
include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF_TRNA      } from '../../modules/local/uncompress_gff'
include { UNCOMPRESS_GFF as UNCOMPRESS_HOST_GFF_TRNA_FILE  } from '../../modules/local/uncompress_gff'

include { COMBINE_FILES as COMBINE_HOST_GENOME_TRNA_GFF_STAR_SALMON;
        COMBINE_FILES as COMBINE_PATHOGEN_HOST_GFF_FILES_HTSEQ; }  from '../../modules/local/combine_files'

include { 
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_HOST_HTSEQ;
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_HOST_SALMON;
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_PATHOGEN_HTSEQ; 
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_PATHOGEN_SALMON 
} from '../../modules/local/extract_annotations'

include { 
    REPLACE_ATTRIBUTE_GFF_STAR_SALMON as REPLACE_ATTRIBUTE_GFF_STAR_SALMON_TRNA;
    REPLACE_ATTRIBUTE_GFF_STAR_SALMON as REPLACE_ATTRIBUTE_GFF_STAR_SALMON_TRNA_FILE;
    REPLACE_ATTRIBUTE_GFF_STAR_SALMON as REPLACE_ATTRIBUTE_GFF_STAR_SALMON_PATHOGEN;
} from '../../modules/local/replace_attribute'


include {
    REPLACE_GENE_FEATURE_GFF_SALMON as REPLACE_GENE_FEATURE_GFF_PATHOGEN_SALMON;
    REPLACE_GENE_FEATURE_GFF_SALMON as REPLACE_GENE_FEATURE_GFF_HOST_SALMON
 } from '../../modules/local/replace_gene_feature'


workflow ANNOTATE_DATA {

    take: 
    
        uncompress_host_gff_trna_file
        uncompress_host_gff_trna
        uncompress_host_gff
        uncompress_pathogen_gff
        
    main:
        // pathogen

        EXTRACT_ANNOTATIONS_PATHOGEN_HTSEQ (
            uncompress_pathogen_gff,
            params.gene_feature_gff_to_quantify_pathogen,
            params.pathogen_gff_attribute,
            params.extract_annotations_pathogen_htseq_organism,
            'htseq'
        )

        REPLACE_ATTRIBUTE_GFF_STAR_SALMON_PATHOGEN(
            uncompress_pathogen_gff,
            params.gene_attribute_gff_to_create_transcriptome_pathogen
        )

        EXTRACT_ANNOTATIONS_PATHOGEN_SALMON (
            REPLACE_ATTRIBUTE_GFF_STAR_SALMON_PATHOGEN.out,
            params.gene_feature_gff_to_create_transcriptome_pathogen,
            "parent",
            params.extract_annotations_pathogen_salmon_organism,
            'salmon'
        )
        // host

        EXTRACT_ANNOTATIONS_HOST_HTSEQ (
            uncompress_host_gff,
            params.gene_feature_gff_to_quantify_host,
            params.host_gff_attribute,
            params.extract_annotations_host_htseq_organism,
            'htseq'
        ) 

        if(params.gff_host_tRNA){
            REPLACE_ATTRIBUTE_GFF_STAR_SALMON_TRNA(
                uncompress_host_gff_trna,
                params.gene_attribute_gff_to_create_transcriptome_host
            )
            REPLACE_ATTRIBUTE_GFF_STAR_SALMON_TRNA_FILE(
                uncompress_host_gff_trna_file,
                params.gene_attribute_gff_to_create_transcriptome_host
            )

            COMBINE_HOST_GENOME_TRNA_GFF_STAR_SALMON(
                REPLACE_ATTRIBUTE_GFF_STAR_SALMON_TRNA.out.mix(
                    REPLACE_ATTRIBUTE_GFF_STAR_SALMON_TRNA_FILE.out
                ).collect()
            )
            REPLACE_GENE_FEATURE_GFF_HOST_SALMON(
                COMBINE_HOST_GENOME_TRNA_GFF_STAR_SALMON.out,
                params.gene_feature_gff_to_create_transcriptome_host
            )

            EXTRACT_ANNOTATIONS_HOST_SALMON (
                REPLACE_GENE_FEATURE_GFF_HOST_SALMON.out,
                params.extract_annotations_host_salmon_feature,
                params.extract_annotations_host_salmon_attribute,
                params.extract_annotations_host_salmon_organism,
                'salmon'
            )   
        }
}

