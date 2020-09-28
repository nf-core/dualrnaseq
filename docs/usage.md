# nf-core/dualrnaseq: Running the pipeline

## Table of contents

1. [Running the pipeline](#1-running-the-pipeline)
   * [Quick start](#11-quick-start)
   * [Basic run](#12-basic-run)
   * [Updating the pipeline](#13-updating-the-pipeline)
   * [Reproducibility](#14-reproducibility)
2. [Configuration profile](#2-configuration-profile)
3. [Input sequence reads](#3-input-sequence-reads)
4. [Reference genomes and annotation](#4-reference-genomes-and-annotation)
   * [Genomes](#41-genomes)
   * [Annotation](#42-annotation)
5. [Read mapping and quantification](#5-read-mapping-and-quantification)
   * [Salmon - Selective alignment](#51-salmon---selective-alignment)
   * [Salmon - quantification in alignment-based mode](#52-salmon---quantification-in-alignment-based-mode)
   * [STAR - alignment-based genome mapping](#53-star---alignment-based-genome-mapping)
6. [Mapping statistics](#6-mapping-statistics)
7. [Example usage](#7-example-usage)
8. [Output files](#8-output-files)
9. [Job resources](#9-job-resources-and-submission)

## 1. Running the pipeline

### 1.1 Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/dualrnaseq -profile test,<docker/singularity/conda/institute>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run nf-core/dualrnaseq -profile <docker/singularity/conda/institute> --reads '*_R{1,2}.fastq.gz' --genome_host GRCh38 --genome_pathogen SL1344
```

See [parameters docs](docs/parameters.md) for all of the available parameters when running the pipeline.

### 1.2 Basic run

Once ready, a basic command for running the pipeline would be the following:

```bash
nextflow run nf-core/dualrnaseq/main.nf -profile docker \
--reads "/folder_to_reads/*_R{1,2}.fq.gz" \
--fasta_host host.fa --fasta_pathogen pathogen.fa \
--gff_host host.gff --gff_pathogen pathogen.gff \
--run_star --outdir results
```

This will launch the pipeline with the `docker` configuration profile (click [here](#https://nf-co.re/usage/configuration), or see [below](#2-configuration-profile) for more information about profiles).
It takes all compressed .fq files in the specified folder and will map the reads (using STAR) to the supplied host and pathogen genomes.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and logs.
```

### 1.3 Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. Subsequent runs will use the cached version if available - even if the pipeline has since been updated. To make sure that you're running the latest version, make sure that you regularly run this command so the cached version is the most recent:

```bash
nextflow pull nf-core/dualrnaseq
```

### 1.4 Reproducibility

It's a good idea to specify a pipeline version when running on your data. This will ensure results can be reproduced more easily, and can be easily referenced or referred to, particularly when multiple versions are available. If a version number was't strictly defined in a previous run, it can be found in the various reports in the `results` directory.

To find the latest version number, go to the [nf-core/dualrnaseq releases page](https://github.com/nf-core/dualrnaseq/releases) - which will be a numeric value (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `nextflow run -r 1.3.1 nf-core/dualrnaseq/main.nf`.

## 2. Configuration profile

**`-profile`**

Use this parameter to choose a configuration profile, which can define specific presets for different compute environments.

Several generic profiles are bundled with the pipeline, instructing it to use software packaged within different container-based objects (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` (the order of arguments is important!). They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

The nf-core/dualrnaseq pipeline contains Docker, Singluarity, Conda and test configuration profiles:

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/dualrnaseq`](http://hub.docker.com/r/nfcore/dualrnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from SingularityHub: [`nfcore/dualrnaseq`](https://singularity-hub.org/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data, requiring no other parameters

## 3. Input sequence reads

`--reads`

Input files can be read as either .fastq or .fastq.gz. They should be named descriptively without spaces and special characters (such as : and @), with the corresponding replicate (if any) denoted with a capital `R` or lower case `r`, and read number `1` or `2` appended at the end. The best practise for this pipeline is to use underscores to separate different experimental conditions, for example:

| Paired-end | Single-end |
| --- | --- |
| Host_pathogen_mock_10h_R1_1.fq | Host_pathogen_mock_10h_r1.fq |
| Host_pathogen_mock_10h_R1_2.fq | Host_pathogen_mock_10h_r2.fq |
| Host_pathogen_mock_10h_R2_1.fq | Host_pathogen_mock_10h_r3.fq |
| Host_pathogen_mock_10h_R2_2.fq | Host_pathogen_MOI_1_10h_R1.fq |
| Host_pathogen_mock_10h_R3_1.fq | Host_pathogen_MOI_1_10h_R2.fq |
| Host_pathogen_mock_10h_R3_2.fq | Host_pathogen_MOI_1_10h_R3.fq |
| Host_pathogen_MOI_1_10h_r1_1.fq |  |
| Host_pathogen_MOI_1_10h_r1_2.fq |  |
| Host_pathogen_MOI_1_10h_r2_1.fq |  |
| Host_pathogen_MOI_1_10h_r2_2.fq |  |
| Host_pathogen_MOI_1_10h_r3_1.fq |  |
| Host_pathogen_MOI_1_10h_r3_2.fq |  |

Once correctly named, instead of typing all files, the folder containing these files can be specified, such as `--reads "folder_to_files/*.fq.gz"`. Or, for paired-end `--reads "folder_to_files/*{1,2}.fastq.gz"`

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### Additional parameters

The following two parameters are associated with the sequencing library type. Their defaults are shown below, but should be changed to the experiment-specific values.

`--single_end`

`--stranded "yes"`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when launched.

> Note: it is not possible to run a mixture of single-end and paired-end files in one run.

## 4. Reference genomes and annotation

### 4.1 Genomes

The main goal of Dual RNA-seq is simultaneous profiling of host and pathogen gene expression. Thus, the pipeline requires references for each of the organisms.

These parameters can be used in two ways:

#### A) Using a configuration file

You can create your own configuration file with sets of reference files and save it here: `...conf/genomes.config`

> See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions.

The syntax for this reference configuration would be the following:

```nextflow
params {
  genomes {
    'GRCh38' {
      fasta_host  = '<path to the genome fasta file>'
      gff_host = '<path to the genome annotation file>'
      gff_host_tRNA = '<path to the tRNA annotation file>' // Optional
      transcriptome_host = '<path to the transcriptome fasta file>' // Optional
             }
    'SL1344' {
      fasta_pathogen  = '<path to the genome fasta file>'
      gff_pathogen = '<path to the genome gff annotation file>'
      transcriptome_pathogen = '<path to the transcriptome fasta file>' // Optional
             }
          }
        // Default genomes (optional). Ignored if --genome_host 'OTHER-GENOME' and --genome_pathogen 'OTHER-GENOME' specified on command line
      genome_host = 'GRCm38'
      genome_pathogen = 'SL1344'
        }
```

Defining default genomes in your configuration file is optional. You can specify the references using the following flags on command line:

`--genome_pathogen SL1344`

`--genome_host GRCm38`

> Any number of additional genomes can be added to this file and specified through either `--genome_host` or `--genome_pathogen`.

If using a custom genome file, you will also need to include either the following line, or something similar in your `nextflow.config` file to make sure the information is being read when the pipeline runs.

```bash
includeConfig 'conf/custom_genomes.config'
```

Note:

* The transcriptome fasta file is created by default in the pipeline using the provided genome and annotation files. If you already have one, you can specify it here as shown above, and through the parameter ```--read_transcriptome_fasta_host_from_file``` or
```--read_transcriptome_fasta_pathogen_from_file```

* If `gff_host_tRNA` file is provided, the pipeline combines `gff_host` and `gff_host_tRNA` files to create host gff file.

* You don't have to specify the path to the host and pathogen transcriptomes in your conf/genomes.config file, as this will be created if needed.

#### B) Using pipeline-specific parameters

If preferred, you can specify each parameter manually and link to appropriate files.

Host:

`--fasta_host` 'path to file'

`--gff_host`

`--gff_host_tRNA`

`--read_transcriptome_fasta_host_from_file`

`--transcriptome_host`

Pathogen:

`--fasta_pathogen`

`--gff_pathogen`

`--read_transcriptome_fasta_pathogen_from_file`

`--transcriptome_pathogen`

> Note: Since many dual RNA-seq experiments are likely to use pathogen-based references that have to be manually downloaded. We recommend adding a new entry to the `genomes.conf` file as depicted [above](#4-reference-genomes-and-annotation), or through specific parameters of `--fasta_pathogen` and `--gff_pathogen`.

##### Host tRNA

We have specified this parameter for users familiar with the [Gencode gene annotations](https://www.gencodegenes.org/). Their annotative files include lncRNAs, snoRNAs, rRNAs and other non coding genes except tRNAs. tRNAs are available in another gff file (predicted tRNA genes). The tRNA gff file looks a little different than the main annotation file, so we don't recommend adding other gff file in the tRNA gff file path place.

##### Warning! The nf-core/dualrnaseq pipeline does not support iGenomes

Many nf-core pipelines provide a possibility to use iGenomes [Reference genomes](https://nf-co.re/usage/reference_genomes). However, the nf-core/dualrnaseq pipeline does not support this functionality.

### 4.2 Annotation

Host-based annotations are generally more uniform in design than pathogen annotations, where different terms for the same feature are often used. For example, in Human annotation files, main features are defined as genes, transcripts and exons, and associated with gene_id, transcript_id and other uniform identifiers. When extracting features, this uniform naming convention makes this straight forward. However, bacterial naming conventions are less uniform. Names for features include genes, CDS, ID, Name, locus_tag amongst others.

We have defined four parameters for both the host and pathogen to account for these differences in conventions - and should be modified when needed.
Default values are shown below:

Host:

`--gene_attribute_gff_to_create_transcriptome_host "transcript_id"`
  
`--gene_feature_gff_to_create_transcriptome_host ["exon", "tRNA"]`
  
`--gene_feature_gff_to_quantify_host ["exon","tRNA"]`
  
`--host_gff_attribute "gene_id"`
  
Pathogen:

`--gene_attribute_gff_to_create_transcriptome_pathogen "locus_tag"`
  
`--gene_feature_gff_to_create_transcriptome_pathogen ["gene","sRNA","tRNA","rRNA"]`
  
`--gene_feature_gff_to_quantify_pathogen ["gene", "sRNA", "tRNA", "rRNA"]`

`--pathogen_gff_attribute "locus_tag"`

#### Features within pathogen annotation files

As previously mentioned, bacterial annotative files can be challenging to work with due to the non-uniformity. We find that summarising the contents of the GFF can help to identify which features you may want to examine. The following bash script can do this easily:

```bash
awk -F '\t' '{print $3}' file.gff3 | sort | uniq -c
```

## 5. Read mapping and quantification

The nf-core/dualrnaseq pipeline provides three strategies to map and quantify your dual RNA-seq data.

1) The first approach **Salmon with Selective alignment**, performs mapping using the Selective alignment algorithm which quantifies the reads. `--run_salmon_selective_alignment`

2) The second approach **Salmon with alignment-based mode**, utalises aligned reads from STAR to quantify input reads. `--run_salmon_alignment_based_mode`

3) The third strategy utilises alignment-based mapping executed with **STAR**, counting uniquely aligned reads with **HTSeq**. `--run_star` and `--run_htseq_uniquely_mapped`

### 5.1 Salmon - selective alignment

Salmon is a transcriptome-based mapping tool that performes both mapping and quantification. In the first phase it performes indexing of reference transcripts (pathogen transcripts are defined as gene or CDS features), where a chimeric transcriptome of host and pathogen files is created. During this step, coordinates of gene features are also extracted from the host and pathogen annotation files.
To avoid spurious mapping of reads that originate from unannotated locus to sequences similar to annotatated transcripts, a decoy-aware transcriptome is created and incorporated into the index. In the pipeline the decoy sequence is created from both host and pathogen genomes (which are both required).

> Note: **Selective alignment** is an improvement to the original alignment-free approach (also called quasi-mapping). In Selective-alignment, the best transcript for a read from a set of mappings is selected based on the alignment-based score instead of the the longest exact match - which increases accuracy. See [`Salmon documentation`](https://salmon.readthedocs.io/en/latest/salmon.html) for more information.

### 5.2 Salmon - quantification in alignment-based mode

In this [mode](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-alignment-based-mode), Salmon performs quantification utilising an aligned BAM file. In the nf-core/dualrnaseq pipeline, the alignment file is generated with STAR. The first step involves creating an index of a chimeric genome (created from the host and pathogen genome fasta files). Next, STAR performs an alignment, but for the purpose of Salmon (it generates alignments translated into transcript coordinates). To learn more on this behavior, please see `Output in transcript coordinates` from the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

> Note: there are numerous STAR-based flags that can be modified within the pipeline - which can be viewed [here](https://github.com/BarquistLab/nf-core-dualrnaseq/blob/master/docs/parameters.md).
> Salmon performs quantification based on a reference transcriptome. It is recommended to allow the pipeline to create a transcriptome using the provided genome (fasta) and annotative (gff) files.
> When quantifying alignments, the parameters `--libtype` and `--incompatPrior` should be adjusted as required.

### 5.3 STAR - alignment-based genome mapping

STAR is a splice-aware alignment tool which aligns reads to a reference genome. In the nf-core/dualrnaseq pipeline, STAR generates a chimeric genome index, then identifies and maps spliced alignments across splice junctions. Therefore, the paths to host and pathogen genomes and annotative files must be provided either through iGenomes, or directly using the appropriate parameters:  `--genome_host`, `--genome_pathogen`, `--gff_host` and `--gff_pathogen`.

## 6. Mapping statistics

To summarise the mapping statistics including total mapped reads, unmapped reads, host-specifc and pathogen-specific mapped reads, add the following parameter to the command line: `--mapping_statistics`.

This will create the following:

* Count the total number of reads before and after trimming
* Scatterplots comparing all replicates (separate for both host and pathogen reads)
* Plots of the % of mapped/quantified reads
* Plots of RNA-class statistics (for more information click [here](https://github.com/BarquistLab/nf-core-dualrnaseq/blob/master/docs/parameters.md#9-rna-mapping-statistics).

## 7. Example usage

As discussed above in the [Read mapping and quantification section](#52-salmon---quantification-in-alignment-based-mode) above, there are different ways to run the pipeline:

### Example 1

* Using Docker

* Single end reads

* Host (Human) and pathogen (*E. coli*) genomes defined through `genomes.conf`

* Salmon - Selective alignment

```bash
nextflow run nf-core-dualrnaseq/main.nf" -profile docker,cluster \
--genome_host "GRCh38" --genome_pathogen "Escherichia_coli_K_12_DH10B" \
--reads "folder_to_reads/*.fq.gz" --single_end \
--outdir "/outdir_folder/" \
--run_salmon_selective_alignment \
```

### Example 2

* Using Docker

* Paired-end reads (unstranded)

* Host (Mouse) and pathogen (*C. trachomatis*) genomes defined through `genomes.conf`

* Salmon - quantification in alignment-based mode

* Custom kmer length

 ```bash
qsub -q all.q nextflow run nf-core-dualrnaseq/main.nf" -profile docker \
--genome_host "GRCm38" --genome_pathogen "C_trachomatis_strain_d" \
--reads "folder_to_reads/*{1,2}.fastq.gz" \
--outdir "/outdir_folder/" \
--run_salmon_alignment_based_mode --libtype "IU" --incompatPrior 0.0 --kmer_length 19 \
```

### Example 3

* Using Singularity

* Single end reads

* Host (Human) and pathogen (*M. pneumoniae*) genomes defined through `genomes.conf`

* STAR - alignment-based genome mapping

 ```bash
nextflow run nf-core-dualrnaseq/main.nf" -profile singularity \
--genome_host "GRCH38" --genome_pathogen "Mycoplasma_pneumoniae" \
--reads "folder_to_reads/*.fq.gz" --single_end \
--outdir /outdir_folder/ \
--run_star --run_htseq_uniquely_mapped \
```

### Example 4

* Run on a SGE cluster

* Using Singularity

* Single end reads

* Host (Human) and pathogen (*S. typhimurium*) genomes defined through `genomes.conf`

* All three modes with custom gene attributes

```bash
qsub -q all.q nextflow run nf-core-dualrnaseq/main.nf" -profile singularity,cluster \
--genome_host "GRCH38" --genome_pathogen "Salmonella_typhimurium" \
--reads "folder_to_reads/*.fq.gz" --single_end \
--outdir "/outdir_folder/" \
--run_salmon_alignment_based_mode --libtype "SF" \
--run_salmon_selective_alignment \
--run_star --run_htseq_uniquely_mapped \
--gene_feature_gff_to_create_transcriptome_pathogen "[ID, sRNA, tRNA, rRNA]" \
--gene_feature_gff_to_quantify_pathogen "[locus_tag, sRNA, tRNA, rRNA]"
```

## 8. Output files

Click [here](https://github.com/BarquistLab/nf-core-dualrnaseq/blob/master/docs/output.md) for a description on output files.

## 9. Job resources and submission

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus, Nextflow processes must run until the pipeline is finished. To achieve this, we recommend running in the background through `screen` / `tmux` or a similar tool. Alternatively, you can run dualrnaseq submitted by your job scheduler on a cluster.

It is also recommended to limit the Nextflow Java virtual machines memory. This can be achieved by adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

### Automatic resubmission

Each step in the pipeline has a default set of requirements for the number of CPUs, memory and time. For most steps within the pipeline, if the job exits with an error code of `143` (exceeded requested resources), it will automatically resubmit with higher requirements (2 x original, then 3 x original). If these requests continue to fail, then the pipeline will stop.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly, it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

### AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

`--awsqueue`  The JobQueue that you intend to use on AWS Batch.

`--awsregion` The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

`--awscli`  The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.
