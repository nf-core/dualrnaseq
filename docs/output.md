# nf-core/dualrnaseq: Output

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [FastQC](#fastqc)
* [Trimming reads](#trimming-reads)
* [Mapping/Quantification](#mappingquantification)
* [MultiQC](#multiqc)
* [Pipeline info](#pipeline-info)
* [Mapping statistics](#mapping-statistics)

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads.
It provides information about the quality score distribution across your reads and the per base sequence content (%T/A/G/C).
Information about adapter contamination and other overrepresented sequences is also displayed.

**Output directory for raw reads:** `results/fastqc`

**Output directory for trimmed reads:** `results/fastqc_after_trimming`

Contents:

* `sample_fastqc.html` and `sample_trimmed_fastqc.html`
  * FastQC reports,`sample_fastqc.html` and  `sample_trimmed_fastqc.html` contain quality metrics for your untrimmed raw and trimmed fastq files, respectively. 
* `zips/sample_fastqc.zip` and `zips/sample_trimmed_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

## Trimming reads

Two software options are provided to remove low quality reads and adapters: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

All samples are processed in the output directory and appended with `_trimmed.fastq.gz`

**Output directory:** `results/trimming`

## Mapping/quantification

Depending on which mapping/quantification mode was selected, will depend on which folders and files appear here.

### A) Salmon selective alignment

**Output directory:** `results/salmon`

Contents:
* `transcripts_index`
  * All files produced by Salmon in the indexing phase.
* subfolders named as samples.
  * All files and folders produced by Salmon in the quantification step of Selective Alignment mode. You can learn more about salmon outputs in [`Salmon documentation`](https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats). Each subfolder also contains separated quantification results from the host and pathogen stored in `host_quant.sf` and `pathogen_quant.sf` files, respectively. 
* `combined_quant.tsv` 
  * Tab delimited file containing combined quantification results of all samples processed by the pipeline. 
* `host_quant_salmon.tsv` and `pathogen_quant_salmon.tsv`
  * Tab delimited files containing combined quantification results for either host or pathogen from all samples.
* `host_combined_quant_annotations.tsv` and `pathogen_combined_quant_annotations.tsv`
  * Tab delimited files containing quantification results for either host or pathogen and annotations extracted from gff files including transcript_id, transcript_name, gene id, gene_name and	gene_type.
* `host_combined_gene_level.tsv`
  * Host gene-level estimates obtained using tximport.
* `host_combined_quant_gene_level_annotations.tsv`
  * Tab delimited file containing host gene-level estimates and annotations extracted from gff files including gene id, gene_name and	gene_type.


### B) STAR + Salmon - alignment based

**Output directory:** `results/STAR_for_salmon`

Contents:
* `index`
  * All files produced by STAR in the indexing phase.
* subfolders named as samples.
  * All files produced by STAR in the alignment step including `sample_Aligned.sample_toTranscriptome.out.bam`. See `Output in transcript coordinates` from the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

**Output directory:** `results/salmon_alignment_mode`

Contents:
* subfolders named as samples
  * All files and folders produced by Salmon in the quantification step of Selective Alignment mode. You can learn more about salmon outputs in [`Salmon documentation`](https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats). Each subfolder also contains separated quantification results from the host and pathogen stored in `host_quant.sf` and `pathogen_quant.sf` files, respectively. 
* `combined_quant.tsv` 
  * Tab delimited file containing combined quantification results of all samples processed by the pipeline. 
* `host_quant_salmon.tsv` and `pathogen_quant_salmon.tsv`
  * Tab delimited files containing combined quantification results for either host or pathogen from all samples.
* `host_combined_quant_annotations.tsv` and `pathogen_combined_quant_annotations.tsv`
  * Tab delimited files containing quantification results for either host or pathogen and annotations extracted from gff files including transcript_id, transcript_name, gene id, gene_name and	gene_type.
* `host_combined_gene_level.tsv`
  * Host gene-level estimates obtained using tximport.
* `host_combined_quant_gene_level_annotations.tsv`
  * Tab delimited file containing host gene-level estimates and annotations extracted from gff files including gene id, gene_name and	gene_type.

### C) STAR + HTSeq

**Output directory:** `results/STAR`

Contents:
* `index`
  * All files produced by STAR in the indexing phase.
* subfolders named as samples.
  * All files produced by STAR in the alignment step.
* `multimapped_reads`.
  * This folder contains both `sample_cross_mapped_reads.txt` file with list of cross-mapped reads between host and pathogen and `sample_no_crossmapped.bam` file with alignment without cross-mapped reads. 

**Output directory:** `results/HTSeq`

Contents:
* `sample_count_u_m.txt`
  * Quantification results for a sample.
* `quantification_results_uniquely_mapped.tsv` 
  * Tab delimited file containing combined quantification results of all samples processed by the pipeline. 
* `quantification_results_uniquely_mapped_NumReads_TPM.tsv` 
  * Tab delimited file containing HTSeq quantification results and TPM values estimated for each gene in each sample.
* `quantification_stats_uniquely_mapped.tsv` 
  * Satistics extracted from HTSeq quantification results.
* `host_quantification_uniquely_mapped_htseq.tsv` and `pathogen_quantification_uniquely_mapped_htseq.tsv`
  * tab delimited files containing combined quantification results for either host or pathogen from all samples.
* `host_combined_quant_annotations.tsv` and `pathogen_combined_quant_annotations.tsv`
  * tab delimited files containing quantification results for either host or pathogen and annotations extracted from gff files including gene_id, gene_name, gene_type and gene length.

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project.

Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory:** `results/MultiQC`

Contents:
* `ProjectName_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `ProjectName_multiqc_report_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline
* `multiqc_plots/`
  * Directory containing plots from within the .html report, saved as `PDF`, `PNG` and `SVG`

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## References



## Mapping statistics

Depending on which mapping/quantification mode was selected, will depend on which folders and files appear here.
In general, files and images within these folders show the number of reads, mapping statistics for both the host and pathogen, RNA-class statistics, as well as gene and transcript-based metrics.

**Output directory:** `results/mapping_statistics`

Contents:

 **STAR:** `results/STAR`
* `star_mapping_stats.tsv`
   * Tab delimited file containg mapping statistics collected from all samples.
* `mapping_stats_samples_total_reads.tsv`
  * Set of mapping and quantification statistics extracted from `star_mapping_stats.tsv` table used to create `mapping_stats_samples_total_reads.pdf` plot. 
* `mapping_stats_samples_total_reads.pdf` 
  * Visualisation of mapping statistics from `mapping_stats_samples_total_reads.tsv` table.

  ![mapping_stats_samples_total_reads](images/mapping_stats_samples_total_reads_star.png)
* `mapping_stats_samples_percentage.tsv`
  * Mapping statistics from `mapping_stats_samples_total_reads.tsv` table expressed in percentage.
* `mapping_stats_samples_percentage.pdf`
  * Visualisation of mapping statistics from `mapping_stats_samples_percentage.tsv` table.

  ![mapping_stats_samples_percentage](images/mapping_stats_samples_percentage_star.png)
  
**HTSeq:** `results/HTSeq`

* `scatter_plots`
  * Scatter plots showing correlations between TPM values of replicates within the same conditions. The pearson correlation coefficient is calculated using untransformed data. 

  ![scatter_plot_pathogen](images/scatter_plot_pathogen_HTseq.png =638x651x20)
  ![scatter_plot_host](images/scatter_plot_host_HTseq.png)
* `htseq_uniquely_mapped_reads_stats.tsv`
   * Colection of mapping and quantification statistics including number of reads uniquely and multi-mapped to either host or pathogen using STAR, number of cross-mapped reads, unmapped reads, trimmed reads and number of assigned reads to the pathogen and host by HTSeq.
* `mapping_stats_samples_total_reads.tsv`
  * Set of mapping and quantification statistics extracted from `htseq_uniquely_mapped_reads_stats.tsv` table used to create `mapping_stats_samples_total_reads.pdf` plot. 
* `mapping_stats_samples_total_reads.pdf`
  * Visualisation of mapping statistics from `mapping_stats_samples_total_reads.tsv` table.

  ![mapping_stats_star_samples_total_reads](images/mapping_stats_samples_total_reads_htseq.png)
* `mapping_stats_samples_percentage.tsv`
  * Mapping and quantification statistics from `mapping_stats_samples_total_reads.tsv` table expressed in percentage.
* `mapping_stats_samples_percentage.pdf`
  * Visualisation of mapping statistics from `mapping_stats_samples_percentage.tsv` table.

  ![mapping_stats_samples_percentage](images/mapping_stats_samples_percentage_htseq.png)
* `RNA_classes_pathogen`
  * `pathogen_RNA_classes_sum_counts_htseq.tsv`
    * Tab delimited file containing pathogen RNA class statistics for each sample. Sum of number of reads assigned to genes which belong to a specific RNA class.
  * `pathogen_RNA_classes_percentage_htseq.tsv`
    * Pathogen RNA class statistics from `pathogen_RNA_classes_sum_counts_htseq.tsv` table expressed in percentage.
  * `RNA_class_stats_combined_pathogen.pdf`
    * Plot showing pathogen RNA class statistics for all samples from `pathogen_RNA_classes_percentage_htseq.tsv` table. 

    ![RNA_class_stats_combined_pathogen](images/RNA_class_stats_combined_pathogen_htseq.png)
  * `sample.pdf`
    * Visualization of pathogen RNA classes statistics for a sample.

    ![RNA_class_stats_sample_pathogen](images/RNA_class_stats_sample_pathogen_htseq.png)
* `RNA_classes_host`
  * `host_RNA_classes_sum_counts_htseq.tsv`
    * Tab delimited file containing host RNA class statistics for each sample. 
  * `host_RNA_classes_percentage_htseq.tsv`
    * Host RNA class statistics from `host_RNA_classes_sum_counts_htseq.tsv` table expressed in percentage.
  * `host_gene_types_groups_gene.tsv`
    * Tab delimited file containing list of host genes including gene ids, gene names, counts obtained from quantification, and gene type assigned to each gene considering RNA class groups defined by `--RNA_classes_to_replace_host`. For more information check [parameters.md](https://github.com/BarquistLab/nf-core-dualrnaseq/blob/master/docs/parameters.md). This table is created only for examination of results. 
  * `RNA_class_stats_combined_host.pdf`
    * Visualization of host RNA class statistics for all samples from `host_RNA_classes_percentage_htseq.tsv` table. 

    ![RNA_class_stats_combined_host](images/RNA_class_stats_combined_host_htseq.png)
  * `sample.pdf`
    * Visualization of host RNA classes statistics for a sample.

    ![RNA_class_stats_sample_host](images/RNA_class_stats_sample_host_htseq.png)

**Salmon:** `results/salmon`
* `scatter_plots`
  * Scatter plots showing correlations between TPM values of replicates within the same conditions. The pearson correlation coefficient is calculated using untransformed data. 

  ![scatter_plot_pathogen](images/scatter_plot_pathogen_salmon.png)
  ![scatter_plot_host](images/scatter_plot_host_salmon.png)
* `salmon_host_pathogen_total_reads.tsv`
   * Tab delimited file containg mapping statistics collected from all samples.
* `mapping_stats_samples_total_reads.tsv`
  * Set of mapping and quantification statistics extracted from `salmon_host_pathogen_total_reads.tsv` table used to create `mapping_stats_samples_total_reads.pdf` plot. 
* `mapping_stats_samples_total_reads.pdf`
  * Visualisation of mapping statistics from `mapping_stats_samples_total_reads.tsv` table.

  ![mapping_stats_star_samples_total_reads](images/mapping_stats_samples_total_reads_salmon.png)
* `mapping_stats_samples_percentage.tsv`
  * Mapping statistics from `mapping_stats_samples_total_reads.tsv` table expressed in percentage.
* `mapping_stats_samples_percentage.pdf`
  * Visualisation of mapping statistics from `mapping_stats_samples_percentage.tsv` table.

  ![mapping_stats_samples_percentage](images/mapping_stats_samples_percentage_salmon.png)
* `RNA_classes_pathogen`
  * `pathogen_RNA_classes_sum_counts_salmon.tsv`
    * Tab delimited file containing pathogen RNA class statistics for each sample.
  * `pathogen_RNA_classes_percentage_salmon.tsv`
    * Pathogen RNA class statistics from `pathogen_RNA_classes_sum_counts_salmon.tsv` table expressed in percentage.
  * `RNA_class_stats_combined_pathogen.pdf`
    * Plot showing pathogen RNA class statistics for all samples from `pathogen_RNA_classes_percentage_htseq.tsv` table. 

    ![RNA_class_stats_combined_pathogen](images/RNA_class_stats_combined_pathogen_salmon.png)
  * `sample.pdf`
    * Visualization of pathogen RNA classes statistics for a sample.

    ![RNA_class_stats_sample_pathogen](images/RNA_class_stats_sample_pathogen_salmon.png)
* `RNA_classes_host`
  * `host_RNA_classes_sum_counts_salmon.tsv`
    * Tab delimited file containing host RNA class statistics for each sample. 
  * `host_RNA_classes_percentage_salmon.tsv`
    * Host RNA class statistics from `host_RNA_classes_sum_counts_salmonhead.tsv` table expressed in percentage.
  * `host_gene_types_groups_transcript.tsv`
    * Tab delimited file containing list of host transcripts including transcript_id, transcript_name, gene ids, gene names, counts obtained from quantification, and gene type assigned to each gene considering RNA class groups defined by `--RNA_classes_to_replace_host`. For more information check [parameters.md](https://github.com/BarquistLab/nf-core-dualrnaseq/blob/master/docs/parameters.md). This table is created only for examination of results. 
  * `RNA_class_stats_combined_host.pdf`
    * Visualization of host RNA class statistics for all samples from `host_RNA_classes_percentage_salmon.tsv` table. 

    ![RNA_class_stats_combined_host](images/RNA_class_stats_combined_host_salmon.png)
  * `sample.pdf`
    * Visualization of host RNA classes statistics for a sample.

    ![RNA_class_stats_sample_host](images/RNA_class_stats_sample_host_salmon.png)

 **Salmon alignment based:** `results/salmon_alignment_based`
 * `scatter_plots`
  * Scatter plots showing correlations between TPM values of replicates within the same conditions. The pearson correlation coefficient is calculated using untransformed data. 

  ![scatter_plot_pathogen](images/scatter_plot_pathogen_salmon_al.png)
  ![scatter_plot_host](images/scatter_plot_host_salmon_al.png)
* `salmon_alignment_host_pathogen_total_reads.tsv`
   * Colection of mapping and quantification statistics including number of reads uniquely and multi-mapped to either host or pathogen transcriptome using STAR, unmapped reads, trimmed reads and number of assigned reads to the pathogen and host by HTSeq.
* `mapping_stats_samples_total_reads.tsv`
  * Set of mapping and quantification statistics extracted from `salmon_alignment_host_pathogen_total_reads.tsv` table used to create `mapping_stats_samples_total_reads.pdf` plot. 
* `mapping_stats_samples_total_reads.pdf`
  * Visualisation of mapping statistics from `mapping_stats_samples_total_reads.tsv` table.

  ![mapping_stats_star_samples_total_reads](images/mapping_stats_samples_total_reads_salmon_al.png)
* `mapping_stats_samples_percentage.tsv`
  * Mapping statistics from `mapping_stats_samples_total_reads.tsv` table expressed in percentage.
* `mapping_stats_samples_percentage.pdf`
  * Visualisation of mapping statistics from `mapping_stats_samples_percentage.tsv` table.

  ![mapping_stats_samples_percentage](images/mapping_stats_samples_percentage_salmon_al.png)
* `RNA_classes_pathogen`
  * `pathogen_RNA_classes_sum_counts_salmon.tsv`
    * Tab delimited file containing pathogen RNA class statistics for each sample.
  * `pathogen_RNA_classes_percentage_salmon.tsv`
    * Pathogen RNA class statistics from `pathogen_RNA_classes_sum_counts_salmon.tsv` table expressed in percentage.
  * `RNA_class_stats_combined_pathogen.pdf`
    * Plot showing pathogen RNA class statistics for all samples from `pathogen_RNA_classes_percentage_htseq.tsv` table. 

    ![RNA_class_stats_combined_pathogen](images/RNA_class_stats_combined_pathogen_salmon_al.png)
  * `sample.pdf`
    * Visualization of pathogen RNA classes statistics for a sample.

    ![RNA_class_stats_sample_pathogen](images/RNA_class_stats_sample_pathogen_salmon_al.png)
* `RNA_classes_host`
  * `host_RNA_classes_sum_counts_salmon.tsv`
    * Tab delimited file containing host RNA class statistics for each sample. 
  * `host_RNA_classes_percentage_salmon.tsv`
    * Host RNA class statistics from `host_RNA_classes_sum_counts_salmon.tsv` table expressed in percentage.
  * `host_gene_types_groups_transcript.tsv`
    * Tab delimited file containing list of host transcripts including transcript_id, transcript_name, gene ids, gene names, counts obtained from quantification, and gene type assigned to each gene considering RNA class groups defined by `--RNA_classes_to_replace_host`. For more information check [parameters.md](https://github.com/BarquistLab/nf-core-dualrnaseq/blob/master/docs/parameters.md). This table is created only for examination of results. 
  * `RNA_class_stats_combined_host.pdf`
    * Visualization of host RNA class statistics for all samples from `host_RNA_classes_percentage_salmon.tsv` table. 

    ![RNA_class_stats_combined_host](images/RNA_class_stats_combined_host_salmon_al.png)
  * `sample.pdf`
    * Visualization of host RNA classes statistics for a sample.

    ![RNA_class_stats_sample_host](images/RNA_class_stats_sample_host_salmon_al.png)


## Pipeline info

Directory containing all of the pipeline-specific reports, timeline and descriptions.

Detailed descriptions of each file can be found on the [Nextflow website](https://www.nextflow.io/docs/latest/tracing.html)

**Output directory:** `results/pipeline_info`

Contents:
* `execution_report.html`
* `execution_timeline.html`
* `execution_trace.txt`
* `pipeline_dag.svg`
* `pipeline_report.html`
* `pipeline_report.txt`
* `results_description.html`
* `software_versions.csv`