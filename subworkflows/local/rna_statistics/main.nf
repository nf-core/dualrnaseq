include {
    RNA_CLASS_STATISTICS_UNIQUELY_MAPPED as RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_HOST;
    RNA_CLASS_STATISTICS_UNIQUELY_MAPPED as RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_PATOGEN;
} from '../../../modules/local/rna_class_statistics_uniquely_mapped/main.nf'

include {
    PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH as PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_HOST;
    PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH as PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_PATOGEN;
} from '../../../modules/local/plot_rna_class_uniquely_mapped_each/main.nf'


include {
    PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED as PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED_HOST;
    PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED as PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED_PATHOGEN;
} from '../../../modules/local/plot_rna_class_uniquely_mapped_combined/main.nf'


include {
    SCATTER_PLOT as SCATTER_PLOT_HOST;
    SCATTER_PLOT as SCATTER_PLOT_PATHOGEN;
} from '../../../modules/local/scatter_plot/main.nf'

/*
    Common statistics used for Star HTSEQ and SALOMON SELECTIVE ALIGNMENT
*/
workflow RNA_STATISTICS {
    take:
        host_quantification
        pathogen_quantification
        scatterplots_host
        scatterplots_pathogen
        annotations_host_htseq
        annotations_pathogen_htseq
        attribute // gff attribute
        tool // salmon, htseq

    main:
        ch_versions = Channel.empty()
        SCATTER_PLOT_HOST(
            host_quantification,
            attribute,
            Channel.value("replicates"), // output from CHECK_REPLICATES
            scatterplots_host,
            "host"
        )
        ch_versions = ch_versions.mix(SCATTER_PLOT_HOST.out.versions.first())
        SCATTER_PLOT_PATHOGEN(
            pathogen_quantification,
            attribute,
            Channel.value("replicates"), // output from CHECK_REPLICATES
            scatterplots_pathogen,
            "pathogen"
        )
        ch_versions = ch_versions.mix(SCATTER_PLOT_PATHOGEN.out.versions.first())

        // // statistics 
        RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_HOST(
            host_quantification,
            attribute,
            annotations_host_htseq,
            Channel.empty(),
            "host",
            "htseq"
        )

        RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_PATOGEN(
            pathogen_quantification,
            attribute,
            annotations_pathogen_htseq,
            Channel.empty(),
            "pathogen",
            "htseq"
        )

        host_input_for_plots = RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_HOST.out.tsv.filter{it.countLines > 0}
        pathogen_input_for_plots = RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_HOST.out.tsv.filter{it.countLines > 0}
        countLines
        // each
        PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_HOST(
            host_input_for_plots,
            tool
        )
        PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_PATOGEN(
            pathogen_input_for_plots,
            tool
        )

        // combined
        PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED_HOST(
            host_input_for_plots,
            tool
        )
        PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED_PATHOGEN(
            pathogen_input_for_plots,
            tool
        )
    emit:
        versions = ch_versions  // channel: [ versions.yml ]
}