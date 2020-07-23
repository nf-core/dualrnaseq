# nf-core/dualrnaseq: Usage

## Table of contents

* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`-profile`](#-profile)
  * [`--reads`](#--reads)
  * [`--single_end`](#--single_end)
* [Reference genomes](#reference-genomes)
  * [`--genome_host` (using iGenomes)](#--genome_host-using-igenomes)
  * [`--genome_pathogen` (using iGenomes)](#--genome_pathogen-using-igenomes)
  * [`--fasta_host`](#--fasta_host)
  * [`--fasta_pathogen`](#--fasta_pathogen)
  * [`--gff_host`](#--gff_host)
  * [`--gff_host_tRNA`](#--gff_host_tRNA)
  * [`--gff_pathogen`](#--gff_pathogen)
  * [`--read_transcriptome_fasta_host_from_file`](#--read_transcriptome_fasta_host_from_file)
  * [`--read_transcriptome_fasta_pathogen_from_file`](#--read_transcriptome_fasta_pathogen_from_file)
  * [`--transcriptome_host`](#--transcriptome_host)
  * [`--transcriptome_pathogen`](#--transcriptome_pathogen)
  * [`--igenomes_ignore`](#--igenomes_ignore)
* [Fastqc](#Fastqc)
    * [`-skipFastqc`](#--skipFastqc)
* [Adapter trimming](#Adapter-trimming)
    * [`--a`](#--a)
    * [`--A`](#--A)
    * [`--quality-cutoff`](#--quality-cutoff)
    * [`--skipTrimming`](#--skipTrimming)
* [Read mapping and quantification](#Read-mapping-and-quantification)
  * [Salmon - Selective-Alignment](#Salmon---Selective-Alignment)
    * [`--run_salmon_selective_alignment`](#--run_salmon_selective_alignment)
    * [`--gene_feature_gff_to_create_transcriptome_host`](#--gene_feature_gff_to_create_transcriptome_host)
    * [`--gene_feature_gff_to_create_transcriptome_pathogen`](#--gene_feature_gff_to_create_transcriptome_pathogen)
    * [`--gene_attribute_gff_to_create_transcriptome_host`](#--gene_attribute_gff_to_create_transcriptome_host)
    * [`--gene_attribute_gff_to_create_transcriptome_pathogen`](#--gene_attribute_gff_to_create_transcriptome_pathogen)
    * [`--kmer_length`](#--kmer_length)
    * [`--libtype`](#--libtype)
    * [`--writeUnmappedNames`](#--writeUnmappedNames)
    * [`--softclipOverhangs`](#--softclipOverhangs)
    * [`--incompatPrior`](#--incompatPrior)
    * [`--dumpEq`](#--dumpEq)
    * [`--writeMappings`](#--writeMappings)
  * [STAR - alignment-based genome mapping](#STAR---alignment-based-genome-mapping)
    * [`--run_star`](#--run_star)
    * [`--outSAMunmapped`](#--outSAMunmapped)
    * [`--outSAMattributes`](#--outSAMattributes)
    * [`--outFilterMultimapNmax`](#--outFilterMultimapNmax)
    * [`--outFilterType`](#--outFilterType)
    * [`--alignSJoverhangMin`](#--alignSJoverhangMin)
    * [`--alignSJDBoverhangMin`](#--alignSJDBoverhangMin)
    * [`--outFilterMismatchNmax`](#--outFilterMismatchNmax)
    * [`--outFilterMismatchNoverReadLmax`](#--outFilterMismatchNoverReadLmax)
    * [`--alignIntronMin`](#--alignIntronMin)
    * [`--alignIntronMax`](#--alignIntronMax)
    * [`--alignMatesGapMax`](#--alignMatesGapMax)
    * [`--outWigType`](#--outWigType)
    * [`--outWigStrand`](#--outWigStrand)
  * [HTSeq - counting of uniquely-mapped reads](#HTSeq---counting-of-uniquely-mapped-reads)
  * [HTSeq - counting of multi-mapped reads](#HTSeq---counting-of-multi-mapped-reads)
  * [Salmon - quantification in alignment-based mode](#Salmon---quantification-in-alignment-based-mode)
    * [`--run_salmon_alignment_based_mode`](#--run_salmon_alignment_based_mode)
    * [`--quantTranscriptomeBan`](#--quantTranscriptomeBan)
* [Maping statistics](#Maping-statistics)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)


## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO nf-core: Document required command line parameters to run the pipeline-->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/dualrnaseq --reads '*_R{1,2}.fastq.gz' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
More information about running the pipeline can be found in the `docs/` directory:

- [Output and how to interpret the results](output.md) - change !!!!!!!!!!!!!!!!!!!!!!!!!!
- [Extra Documentation on annotation](annotation.md) - change !!!!!!!!!!!!!!!!!!!!

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/dualrnaseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/dualrnaseq releases page](https://github.com/nf-core/dualrnaseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/dualrnaseq`](http://hub.docker.com/r/nfcore/dualrnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/dualrnaseq`](http://hub.docker.com/r/nfcore/dualrnaseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

<!-- TODO nf-core: Document required command line parameters -->

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

## Reference genomes
The main goal of Dual RNA-seq is simultaneous profiling of host and pathogen gene expression. Thus, the pipeline requires to provide references for each of the organisms ([genome_host](#--genome_host-(using-iGenomes)) and [genome_pathogen](#--genome_pathogen-(using-iGenomes))).

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.
  
### `--genome_host` (using iGenomes)
There are 3 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome_host` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common host genomes that are supported are:

* Human
  * `--genome_host GRCh38`
* Mouse
  * `--genome_host GRCm38`

<!--  > There are numerous others - check the config file for more.-->

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh38' {
      fasta_host  = '<path to the genome fasta file>'
      gff_host = '<path to the genome gff annotation file>'
      gff_host_tRNA = '<path to the tRNA gff annotation fasta file>'
      transcriptome_host = '<path to the tRNA gff annotation fasta file>'
    }
  }
// Any number of additional genomes, key is used with --genome_host
}
```
The transcriptome fasta file is created by default in the pipeline using provided fasta files and gff file. You don't have to specify the path to the pathogen transcriptome in your conf/genomes.config file. 
If you want to use transcriptome defined in [iGenomes config file](../conf/igenomes.config) or
your conf/genomes.config file, please specify the [`--read_transcriptome_fasta_host_from_file`](###--read_transcriptome_fasta_host_from_file`) flag. 

If `gff_host_tRNA` file is provided, the pipeline combines `gff_host` and `gff_host_tRNA` files to create host gff file.

### `--genome_pathogen` (using iGenomes)
To run the pipeline with pathogen references defined in iGenomes, you must specify which iGenomes keys to use with the `--genome_pathogen` flag.

You can fing the iGenomes keys in [iGenomes config file](../conf/igenomes.config). Common pahogen genome that is supported is:

* _Salmonella_ Typhimurium
  * `--genome_pathogen SL1344`

If your genome of interest is not provided with iGenomes you can create your own configuration file and save it as conf/genomes.config.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'SL1344' {
      fasta_pathogen  = '<path to the genome fasta file>'
      gff_pathogen = '<path to the genome gff annotation file>'
      transcriptome_pathogen = '<path to the genome gff annotation file>'
    }
  }
// Any number of additional genomes, key is used with --genome_pathogen
}
```
The transcriptome fasta file is created by default in the pipeline using provided genomic fasta and gff file. You don't have to specify the path to the pathogen transcriptome in your conf/genomes.config file.

If you want to use transcriptome defined in [iGenomes config file](../conf/igenomes.config) or
your conf/genomes.config file, please specify the [`--read_transcriptome_fasta_pathogen_from_file`](###--read_transcriptome_fasta_pathogen_from_file`) flag. 

### `--fasta_host`

If you prefer, you can specify the full path to your host genome fasta file when you run the pipeline:

```bash
--fasta_host '[path to genome fasta reference of host]'
```

### `--fasta_pathogen`

If you prefer, you can specify the full path to your pathogen genome fasta file when you run the pipeline:

```bash
--fasta_pathogen '[path to genome fasta reference of pathogen]'
```

### `--gff_host`

If you prefer, you can specify the full path to genome annotation file of your host in the gff3 format:

```bash
--gff_host '[path to host genome gff file]'
```
### `--gff_host_tRNA`

Using tRNA annotation is not required, however it is recommended. 
If you prefer, you can specify the full path to tRNA annotations of your host in the gff3 format:

```bash
--gff_host_tRNA '[path to host tRNA gff file]'
```

### `--gff_pathogen`

If you prefer, you can specify the full path to genome annotations of your pathogen in the gff3 format:

```bash
--gff_pathogen '[path to pathogen genome gff file]'
```

### `--read_transcriptome_fasta_host_from_file`
If you prefer to use host transcriptome file defined in either [`--genome_host` (using iGenomes)](#--genome_host-using-igenomes) or [`--transcriptome_host`](#--transcriptome_host), please, specify the `--read_transcriptome_fasta_host_from_file` flag.

```bash
--read_transcriptome_fasta_host_from_file
```
By default, this flag is set to false and the pipeline creates a transcriptome based on provided fasta files and gff files. 

### `--read_transcriptome_fasta_pathogen_from_file`
If you prefer to use pathogen transcriptome file defined in either [`--genome_pathogen` (using iGenomes)](#--genome_pathogen-using-igenomes) or [`--transcriptome_pathogen`](#--transcriptome_pathogen), please, specify the `--read_transcriptome_fasta_pathogen_from_file` flag.

```bash
--read_transcriptome_fasta_pathogen_from_file
```

By default, this flag is set to false and the pipeline creates a transcriptome based on provided fasta files and gff files. 

### `--transcriptome_host`
If you specify the full path to your host transcriptome fasta file with `--transcriptome_host` and set [`--read_transcriptome_fasta_host_from_file`](###--read_transcriptome_fasta_host_from_file`) flag
the pipeline reads the transcriptome from the file. 

```bash
--transcriptome_host '[path to host ranscriptome fasta file]'
```

### `--transcriptome_pathogen`
If you specify the full path to your host transcriptome fasta file with `--transcriptome_pathogen` and set [`--read_transcriptome_fasta_pathogen_from_file`](###--read_transcriptome_fasta_pathogen_from_file`) flag
the pipeline reads the transcriptome from the file. 

```bash
--transcriptome_pathogen '[path to transcriptome fasta file of pathogen]'
```

### `--igenomes_ignore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.


## Fastqc

### `--skipFastqc`
If you don't want to run Fastqc, please specify the following flag:

```bash
--skipFastqc
```
or set the parameter to true in your config file. 


## Adapter trimming
To remove adapter sequences that were introduced during the library preparation the pipeline utilizes cutadapt.
To learn more on cutadapt and its parameters check the [`cutadapt documentation.`](https://cutadapt.readthedocs.io/en/stable/guide.html) 

By default, the pipeline trims Illumina TruSeq adapters. See [`Illumina TruSeq.`](https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq) 

### `--a`
For single-end reads as well as the first reads of paired-end data, adapter sequence can be specified with `--a` flag. See [`adapter-types.`](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types) 

```bash
--a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
```
### `--A`
For paired-end data, the adapter sequence for the second reads can be defined by `--A`. See [`trimming paired-end reads.`](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads)

```bash
--A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
### `--quality-cutoff`
Cutadapt can also remove low-quality ends of the reads. By default, the 3’ end of each read is trimmed using cutoff 10. You can change this threshold defining the `--quality-cutoff` flag. See [`quality trimming.`](https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming)

```bash
--quality-cutoff 10
```
If you specify two comma-separated cutoffs, the first value represents 5’ cutoff, and the second one 3’ cutoff.
```bash
--quality-cutoff 10,15
```
### `--skipTrimming`
If you don't want to trim your reads, please specify the following flag:
```bash
--skipTrimming
```
or set the parameter to true in your config file. 

# Read mapping and quantification

The nf-core/dualrnaseq provides several strategies to map and quantify your Dual RNA-seq data. The first approach involves [`Salmon with Selective alignment`](#Salmon---Selective-Alignment) which performs mapping using Selective alignment algorith and quantifies the reads. The second strategy utilizes alignment-based mapping executed with [`STAR`](#STAR---alignment-based-genome-mapping). To count aligned reads, the pipeline uses HTSeq in two modes, counting of [`uniquely-mapped reads`](#HTSeq---counting-of-uniquely-mapped-reads) and [`multi-mapped reads`](#HTSeq---counting-of-multi-mapped-reads). The last approach involves [`Salmon with alignment-based mode`](#Salmon---quantification-in-alignment-based-mode). 

## Salmon - Selective-Alignment

Salmon is a transcriptome-based mapping tool that performes both mapping and quantification. In the first phase it performes indexing of reference transcripts. The nf-core/dualrnaseq pipeline creates a chimeric transcriptome of host and pathogen using fasta and gff files defined with [`--genome_host`](#--genome_host-(using-iGenomes)) and [`--genome_pathogen`](#--genome_pathogen-(using-iGenomes)). You can also specify a pathway to each file separatly using [`--fasta_host`](#--fasta_host),[`--fasta_pathogen`](#--fasta_pathogen), [`--gff_host`](#--gff_host), [`--gff_host_tRNA`](#--gff_host_tRNA), and [`--gff_pathogen`](#--gff_pathogen). The pipeline extracts coordinates of gene feautures, defined with [`--gene_feature_gff_to_create_transcriptome_host`](#--gene_feature_gff_to_create_transcriptome_host) and  [`gene_feature_gff_to_create_transcriptome_pathogen`](#gene_feature_gff_to_create_transcriptome_pathogen) from the host and pathogen gff annotation files, respectively.     
If you prefer to use your own transcriptomes, please define the following flags: [`--read_transcriptome_fasta_host_from_file`](#--read_transcriptome_fasta_host_from_file), [`--read_transcriptome_fasta_pathogen_from_file`](#--read_transcriptome_fasta_pathogen_from_file), [`--transcriptome_host`](#--transcriptome_host), and [`--transcriptome_pathogen`](#--transcriptome_pathogen). 
Since Salmon is a transctiptome-based mapping method, to avoid spurious mapping of reads that origanate from unannotated locus with sequece similar to an annotatated transcripts, a decoy-aware transcriptome is created and incorporated into the index. In the pipeline the decoy sequence is created from both pathogen and host entire genomes. Thus, it is required to provide 
host and pathogen genome fasta files either specifying [`--genome_host`](#--genome_host-(using-iGenomes)) and [`--genome_pathogen`](#--genome_pathogen-(using-iGenomes)) or [`--fasta_host`](#--fasta_host) and [`--fasta_pathogen`](#--fasta_pathogen).

Salmon has introduced an improvement to the alignment-free mapping approach, Selective-Alignment. In contrast to quasi-mapping, in selective-alignment mode the best transcript for a read from a set of mappings is selected based on the alignment-based score instead of the the longest exact match which increases the accuracy of the tool. See [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html) 


### `--run_salmon_selective_alignment`

To run Salmon with Selective alignment, please specify the following flag:

```bash
--run_salmon_selective_alignment
```

### `--gene_feature_gff_to_create_transcriptome_host` 
The pipeline uses gene features from the 3rd column of the host gff file to extract the coordinates of transcripts to be quantified. To specify them, please set up the `--gene_feature_gff_to_create_transcriptome_host` flag.

```bash
--gene_feature_gff_to_create_transcriptome_host exon tRNA
```
By default, the pipeline uses `exon` from the [`--gff_host`](#--gff_host) file and `tRNA` from the [`--gff_host_tRNA`](#--gff_host_tRNA) file. 

### `--gene_feature_gff_to_create_transcriptome_pathogen` 
The pipeline uses gene features from the 3rd column of the pathogen gff file to extract the coordinates of transcripts to be quantified. To specify them, please set up the `--gene_feature_gff_to_create_transcriptome_pathogen` flag.
```bash
--gene_feature_gff_to_create_transcriptome_pathogen 
```
By default, the pipeline uses features as `gene`, `sRNA`, `tRNA` and `rRNA` from the [`--gff_pathogen`](#--gff_pathogen) file. 

### `--gene_attribute_gff_to_create_transcriptome_host`
This flag defines a gene attribute of the 9th column of the host gff file which is used to extract the transcript names. 

```bash
--gene_attribute_gff_to_create_transcriptome_host transcript_id
```
By default, the pipeline extracts `transcript_id` from the [`--gff_host`](#--gff_host) file. 

### `gene_attribute_gff_to_create_transcriptome_pathogen`
This flag defines a gene attribute of the 9th column of the pathogen gff file which is used to extract the transcript names. 

```bash
--gene_attribute_gff_to_create_transcriptome_pathogen locus_tag
```
By default, the pipeline extracts `locus_tag` from the [`--gff_pathogen`](#--gff_pathogen) file. 

### `--kmer_length`

To define the k-mer length (`-k` parameter in Salmon, see [`preparing transcriptome indices`](https://salmon.readthedocs.io/en/latest/salmon.html?highlight=index#preparing-transcriptome-indices-mapping-based-mode)) set the `--kmer_length` flag. 

```bash
--kmer_length 21
```
By default, the size of k-mers in the nf-core/dualrnaseq pipeline is set up to 21. 

### `--libtype`

To define the type of sequencing library of your data, specify the following flag:

```bash
--libtype SF
```
To learn more on library types available in Salmon, please read [_`What’s this LIBTYPE?`_](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype)

### `-- writeUnmappedNames`

By default the pipeline saves names of unmapped reads. You can learn more about this option in [`Salmon documentation`](https://salmon.readthedocs.io/en/latest/salmon.html#writeunmappednames). If you don't want to keep this option, set the `--writeUnmappedNames` flag to false.

```bash
--writeUnmappedNames true
```

### `--softclipOverhangs`

By default, the pipeline allows for soft-clipping. 

```bash
--softclipOverhangs true
```
_"Allow soft-clipping of reads that overhang the beginning or ends of the transcript. In this case, the overhaning section of the read will simply be unaligned, and will not contribute or detract from the alignment score"_. 
If it is set to `false`, the end-to-end alignment of the entire read is forced, so that the occurance of overhanings may affect the alignment score.

### `--incompatPrior`

By default, the nf-core/dualrnaseq pipeline set the `--incompatPrior` to 0.0, to ensure that only mappings or alignments that are compatible with the library type are considered by Salmon. You can find more information on this parameter in [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#incompatprior) 

```bash
--incompatPrior
```
### `--dumpEq`
By default, to save the equivalence classes and their counts this option is set to `true`. See [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#dumpeq). 

```bash
--dumpEq true
```
### `--writeMappings`

If set to true, the pipeline will create a `mapping.sam` file with mapping information. To learn more on this option, please look at the [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#writemappings) 

```bash
--writeMappings false
```

## STAR - alignment-based genome mapping

STAR is a splice-aware alignment tool which aligns reads to a reference genome. In the nf-core/dualrnaseq pipeline, STAR generates an index of a chimeric genome. Therefore, the paths to host and pathogen genomes must be provided either through [`--genome_host`](#--genome_host-(using-iGenomes)) and [`--genome_pathogen`](#--genome_pathogen-(using-iGenomes)) or directly using [`--fasta_host`](#--fasta_host) and [`--fasta_pathogen`](#--fasta_pathogen). To identify and correctly map spliced alignments across splice junctions, in the nf-core/dualrnaseq pipeline, STAR utilizes gene annotation provided in gff format using [`--genome_host`](#--genome_host-(using-iGenomes)) and [`--genome_pathogen`](#--genome_pathogen-(using-iGenomes)) or [`--gff_host`](#--gff_host), [`--gff_host_tRNA`](#--gff_host_tRNA), and [`--gff_pathogen`](#--gff_pathogen). 

It genenerates BAM file sorted By coordinates. 

### `--run_star`
To run STAR in the nf-core/dualrnaseq pipeline, please specify the `--run_star` flag.

```bash
--run_star
```
### `--outSAMunmapped`

By default, the pipeline saves unmapped reads within the main BAM file. If you want to switch off this option, set the `--outSAMunmapped` flag to `None`. 

See `--outSAMunmapped` parameter in [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

```bash
--outSAMunmapped Within
```
To save unmapped reads within the main BAM file, please set the `--outSAMunmapped` flag to `Within`. 
By default, this flag is set to `None` in the pipeline. For paired-end reads. or? test...
paired-reads 2nd word: KeepPairs
record unmapped mate for each alignment, and, in case of unsorted
output, keep it adjacent to its mapped mate. Only affects
multi-mapping reads.

### `--outSAMattributes`
To specify the SAM attributes, please use the `--outSAMattributes` flag. To check the possible parameters, you can read more on the `--outSAMattributes` flag in [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

By default, the pipeline uses `Standard` option to keep NH HI AS nM SAM attributes. 

```bash
--outSAMattributes Standard
```

### `--outFilterMultimapNmax`
To specify the maximum number of loci the reads is allowed to map to, please use the following flag:

```bash
--outFilterMultimapNmax 20
```
By default, the option is set to 20 in the pipeline (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)   

### `--outFilterType`

By default, the pipeline keeps reads that contain junctions that passed filtering into
SJ.out.tab. This option reduces the number of ”spurious” junctions. (ENCODE standard options for long RNA-seq pipeline) 

```bash
--outFilterType BySJout
```
You can read more about the flag and its options in the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

### `--alignSJoverhangMin`
To change a number of minimum overhang for unannotated junctions, please set the `--alignSJoverhangMin` flag to a desired number. By default, the pipeline uses 8. (ENCODE standard options for long RNA-seq pipeline) 

```bash
--alignSJoverhangMin 8
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 
### `--alignSJDBoverhangMin`

To set a number of minimum overhang for annotated junctions different than 1 (default option), please specify the --alignSJDBoverhangMin flag. 

```bash
--alignSJDBoverhangMin 1
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

### `--outFilterMismatchNmax`

To define a threshold for a number of mismatches to be allowed, set the `--outFilterMismatchNmax` flag. By default, the pipeline uses a large number to switch off this filter. (ENCODE standard options for long RNA-seq pipeline) 

```bash
--outFilterMismatchNmax 999
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

### `--outFilterMismatchNoverReadLmax`

Using `--outFilterMismatchNoverReadLmax` flag you can define a threshold for a ratio of mismatches to *read* length. The alignment will be considered if the ratio is less than or equal to this value. By default, this option is set to 0.04 in the pipeline. For 2x100b, max number of mismatches is 0.04*200=8 for paired-end reads. (ENCODE standard options for long RNA-seq pipeline)

```bash
--outFilterMismatchNoverReadLmax 0.04
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

### `--alignIntronMin`

By default, the nf-core dualrnaseq pipeline uses 20 as a minimum intron length. If genomic gap is smaller than this value, it is considered as a deletion. 
(ENCODE standard options for long RNA-seq pipeline)

```bash
--alignIntronMin 20
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

### `--alignIntronMax`

The maximum intron length is set to 1000000 in the pipeline (ENCODE standard options for long RNA-seq pipeline). To change this value, please define the desired value for the following flag:

```bash
--alignIntronMax 1000000
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

### `--alignMatesGapMax`
To set the maximum genomic distance between mates, please use the `--alignMatesGapMax` flag. By default, this option is set to 1000000 in the pipeline. (ENCODE standard options for long RNA-seq pipeline)

```bash
--alignMatesGapMax 1000000
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

### `--outWigType`
To generate one of the signal outputs, e.g. "wiggle", "bedGraph", please specify the type of the signal output with the `--outWigType` flag. To find all avaiable types of the singl output, please check the [`STAR documentation`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf). 

By default, the pipeline does not generate any of these files.

```bash
--outWigType None
```
### `--outWigStrand`
To define strandedness of wiggle/bedGraph output, please set the `--outWigStrand` to `Stranded` (default) of `Unstranded`. 

```bash
--outWigStrand Stranded
```
See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

## HTSeq - counting of uniquely-mapped reads

run_star = true

### `--run_htseq_uniquely_mapped`


  gene_feature_gff_to_quantify_host = ["exon","tRNA"]
  gene_feature_gff_to_quantify_pathogen = ["gene", "sRNA", "tRNA", "rRNA"]
  host_gff_atribute = "gene_id"
  pathogen_gff_atribute = "locus_tag"

  tRNA gff - assume it has the same gene atribute as in genomic gff - gene_id
In the nf-core/dualrnaseq pipeline you can specify the following cutadapt parameters: 

## HTSeq - counting of multi-mapped reads

### `--run_htseq_multi_mapped`

  gene_feature_gff_to_quantify_host = ["exon","tRNA"]
  gene_feature_gff_to_quantify_pathogen = ["gene", "sRNA", "tRNA", "rRNA"]
  host_gff_atribute = "gene_id"
  pathogen_gff_atribute = "locus_tag"




## Salmon - quantification in alignment-based mode

In [`alignment-based mode`](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-alignment-based-mode), Salmon performes quantification utilizing an alignment BAM file. In the nf-core/dualrnaseq pipeline the alignment file is generated with STAR. The first step involves indexing of a chimeric genome, creaded from the host and pathogen genome fasta files. 
In the next step, STAR performs an alignment, but on purpose of the Salmon, it generates alignments translated into transcript coordinates. To learn more on this behavior, please check _`Output in transcript coordinates`_ in [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

To run Salmon in alignment-based mode you must specify the following reference files: [`--fasta_host`](#--fasta_host) and [`--fasta_pathogen`](#--fasta_pathogen) as well as [`--gff_host`](#--gff_host), [`--gff_host_tRNA`](#--gff_host_tRNA) (optional), and [`--gff_pathogen`](#--gff_pathogen). You can also define those files in the the config file and then use [`--genome_host`](#--genome_host-(using-iGenomes)) and [`--genome_pathogen`](#--genome_pathogen-(using-iGenomes)) flags. 

You can change the default STAR behaviour defined in the pipeline through the following STAR flags: [`--outSAMunmapped`](#--outSAMunmapped), [`--outSAMattributes`](#--outSAMattributes), [`--outFilterMultimapNmax`](#--outFilterMultimapNmax), [`--outFilterType`](#--outFilterType), [`--alignSJoverhangMin`](#--alignSJoverhangMin), [`--alignSJDBoverhangMin`](#--alignSJDBoverhangMin), [`--outFilterMismatchNmax`](#--outFilterMismatchNmax), [`--outFilterMismatchNoverReadLmax`](#--outFilterMismatchNoverReadLmax), [`--alignIntronMin`](#--alignIntronMin), [`--alignIntronMax`](#--alignIntronMax), and [`--alignMatesGapMax`](#--alignMatesGapMax).  

It genenerates unsorted BAM file. 

Salmon performs quantification for a reference transcriptome and it is recommended to allow the pipeline to create a transcriptome using the provided genome fasta files and annotation gff files. Keep [`--read_transcriptome_fasta_host_from_file`](#--read_transcriptome_fasta_host_from_file) and [`--read_transcriptome_fasta_pathogen_from_file`](#--read_transcriptome_fasta_pathogen_from_file) flags set to false. 

To create the chimeric gff file used by STAR, and the chimeric transcriptome used by Salmon, please define the following flags: [`--gene_feature_gff_to_create_transcriptome_host`](#--gene_feature_gff_to_create_transcriptome_host), [`--gene_feature_gff_to_create_transcriptome_pathogen`](#--gene_feature_gff_to_create_transcriptome_pathogen), [`--gene_attribute_gff_to_create_transcriptome_host`](#--gene_attribute_gff_to_create_transcriptome_host), [`--gene_attribute_gff_to_create_transcriptome_pathogen`](#--gene_attribute_gff_to_create_transcriptome_pathogen)    


To quantify the alignments, please specify the [`--libtype`](#--libtype) and [`--incompatPrior`](#--incompatPrior) Salmon flags. 

### `--run_salmon_alignment_based_mode`
To run Salmon in alignment-based mode, please specify the following flag:

```bash
--run_salmon_alignment_based_mode
```
### `--quantTranscriptomeBan`
The nf-core/dualrnaseq pipeline runs STAR to generate a transcriptomic alignments allowing, by default, for insertions, deletions and soft-clips (`Singleend` option). To prohibit this behaviour, please specify `IndelSoftclipSingleend` for `--quantTranscriptomeBan` flag. See _`--quantTranscriptomeBan`_ in [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 
 
```bash
--quantTranscriptomeBan Singleend
```

## Maping statistics
### `--mapping_statistics`  
def. true

calc. count_total_reads , count_total_trimmed_reads, scatterplots

salmon - extract_processed_reads , unmapped

plot_salmon_mapping_stats_host_pathogen,
RNA class statistics - host, RNA_classes_to_replace.csv , def "data/RNA_classes_to_replace.csv" host
combined, each sample



## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
