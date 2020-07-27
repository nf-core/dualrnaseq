# nf-core/dualrnaseq: Usage

**Table of contents**

1. [Introduction](#introduction)
2. [Running the pipeline](#running-the-pipeline)
   * [Updating the pipeline](#updating-the-pipeline)
   * [Reproducibility](#reproducibility)
3. [Configuration profile](#configuration-profile)   
4. [Input sequence reads](#input-sequence-reads)
5. [Reference genomes and annotation](#reference-genomes-and-annotation)
   * [Genomes](#genomes)
   * [Annotation](#annotation)
6. [Example usage](#example-usage)
7. [Read mapping and quantification](Read-mapping-and-quantification)
8. [Output files](#output-files)
9. [Job resources](#job-resources)
10. [Mapping statistics](#Mapping-statistics)



## 1. Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
 

## 2. Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/dualrnaseq/main.nf -profile docker \
--reads '/folder_to_reads/*_R{1,2}.fq.gz \
--genome_host host.fa.gz --genome_pathogen pathogen.fa.gz \
--run_star --outdir results' 
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.
It takes all compressed .fq files in the specified folder and will map the reads (using STAR) to the supplied host and pathogen genomes. 

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### 2.1 Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. Subsequent runs will use the cached version if available - even if the pipeline has since been updated. To make sure that you're running the latest version, make sure that you regularly run this command so the caches version is the most recent:

```bash
nextflow pull nf-core/dualrnaseq
```

### 2.2 Reproducibility

It's a good idea to specify a pipeline version when running on your data. This will ensure results can be reproduced more easily, and can be easily referenced or referred to, particularly when multiple versions are available. If a version number was't strictly defined in a previous run, it can be found in the various reports in the results directory.

To find the latest version number, go to the [nf-core/dualrnaseq releases page](https://github.com/nf-core/dualrnaseq/releases) - which will be numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.


## 3. Configuration profile

```-profile```

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



## 4. Input sequence reads

```--reads```

Input files can be read as either .fastq or .fastq.gz. They should be named descriptively without spaces and special characters (such as : and @), with the corresponding replicate (if any) appended at the end. The best practise for this pipeline is to use underscores to separate different experimental conditions, for example:

 - Host_pathogen_mock_10h_r1.fq 
 - Host_pathogen_mock_10h_r2.fq 
 - Host_pathogen_mock_10h_r3.fq 
 - Host_pathogen_MOI_1_10h_r1.fq
 - Host_pathogen_MOI_1_10h_r2.fq
 - Host_pathogen_MOI_10_10h_r1.fq
 - Host_pathogen_MOI_10_10h_r2.fq

Once named, instead of typing all files, the folder containing all files can be specified, such as ```--reads "folder_to_files/*.fq.gz"```

Or, for paired-end ```--reads "folder_to_files/*{1,2}.fastq.gz"```


Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`


**Additional parameters**

The following two parameters are associated with the sequencing library type. Their defaults are shown below, but should be changed to the experiment-specific values

```single_end = false```

```stranded = "yes"```
 
> Note: By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when launched. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.


## 5. Reference genomes and annotation

### 5.1 Genomes

The main goal of Dual RNA-seq is simultaneous profiling of host and pathogen gene expression. Thus, the pipeline requires to provide references for each of the organisms

```--genome_host``` 

```--genome_pathogen```

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.


These parameters can be used in three ways: 


**A) Using iGenomes**

There are a range of different species supported with iGenomes references. To run the pipeline, you must specify which one you would like to you with the `--genome_host` flag.

You can find the keys to specify the genomes in the iGenomes config file ```../conf/igenomes.config```. 

Common host genomes that are supported are:

Human: ```--genome_host GRCh38```

Mouse: ```--genome_host GRCm38```

> Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh38' {
      fasta_host  = '<path to the genome fasta file>'
      gff_host = '<path to the genome annotation file>'
      gff_host_tRNA = '<path to the tRNA annotation file>' ## Optional
      transcriptome_host = '<path to the transcriptome fasta file>' ## Optional
            }
          }
        }
```
> Any number of additional genomes can be added to this file and specified through either ```--genome_host``` or ```--genome_pathogen```. 

Note:
 - The transcriptome fasta file is created by default in the pipeline using the provided genome and annotation files. If you already have one, you can specify it here as shown above, or through the parameter ```--read_transcriptome_fasta_host_from_file```
 - If `gff_host_tRNA` file is provided, the pipeline combines `gff_host` and `gff_host_tRNA` files to a create host gff file.
 - You don't have to specify the path to the pathogen transcriptome in your conf/genomes.config file.


**B) Additional configuration file**

If your genome of interest is not provided with iGenomes you can create your own configuration file and save it here: ```...conf/genomes.config```

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'SL1344' {
      fasta_pathogen  = '<path to the genome fasta file>'
      gff_pathogen = '<path to the genome gff annotation file>'
      transcriptome_pathogen = '<path to the transcriptome fasta file>' ## Optional
            }
          }
        }
```

Then to use this reference, within the user-defined set of parameters: ```--genome_pathogen SL1344```  

Note:
 - The transcriptome fasta file is created by default in the pipeline using the provided genome and annotation files. If you already have one, you can specify it here as shown above, or through the parameter ```--read_transcriptome_fasta_pathogen_from_file```


**C) Specific parameters**

If preferred, you can specify each parameter manually and link to appropriate files.


```--fasta_host``` '[path to file]'

```--fasta_pathogen``` 

```--gff_host``` 

```--gff_host_tRNA``` 

```--gff_pathogen``` 

```--read_transcriptome_fasta_host_from_file```

```--read_transcriptome_fasta_pathogen_from_file```

```--transcriptome_host```

```--transcriptome_pathogen```

> Note: If using multiple manual parameters, it may be a good idea to choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

```--igenomes_ignore```




### 5.2 Annotation

Host-based annotations are generally more uniform in design than pathogen annotations, where different terms are often substituted. For example, in Human annotations (.gff or .gtf), main features are defined as genes and exons, and associated with gene_id and transcript_id's. When extracting features, this uniform naming structure makes this straight forward. Although bacterial gff files do not contain exons or transcripts, their naming convensions are less uniform. Names include genes, CDS, ID, Name, locus_tag and many others. 

We have defined four parameters for both the Host and Pathogen to account for these different naming convensions, which can be modified if needed. The defaults are shown below:

Host 

```--gene_attribute_gff_to_create_transcriptome_host = "transcript_id"```
  
```--gene_feature_gff_to_create_transcriptome_host = ["exon", "tRNA"]```
  
```--gene_feature_gff_to_quantify_host = ["exon","tRNA"]```
  
```--host_gff_attribute = "gene_id"```
  
Pathogen

```--gene_attribute_gff_to_create_transcriptome_pathogen = "locus_tag"```
  
```--gene_feature_gff_to_create_transcriptome_pathogen = ["gene","sRNA","tRNA","rRNA"]```
  
```--gene_feature_gff_to_quantify_pathogen = ["gene", "sRNA", "tRNA", "rRNA"]```  

```--pathogen_gff_attribute = "locus_tag"```


We recommend using iGenomes (as discussed above) and the corosponding annotation files (.gff and .gtf) where available. 

However, we realise that many dual RNA-seq experiments are likely to use references that have to be manually downloaded. In these instances, there is no specific parameter to select a specific annotative file and the ```genomes.conf``` file should be edited directly as discussed above. 



## 6. Example usage

Set the default locations
```
dualrnaseq_main="/vol/projects/bmikagos/Dual_RNA_seq/nf-core-dualrnaseq/main.nf"
host_genome="GRCh38"
pathogen_genome="SL1344"
out_dir_dualrnaseq="/vol/projects/bmikagos/Dual_RNA_seq/Hela_Salmonella_nfcore_dual_rna_seq"
libtype="SF"
RNA_classes_to_replace_host="/vol/projects/bmikagos/Dual_RNA_seq/nf-core-dualrnaseq/data/RNA_classes_to_replace.csv"
```
Code to run
```
qsub -q all.q ~/source/nextflow run ${dualrnaseq_main} -profile singularity,cluster \
--genome_host ${host_genome} --genome_pathogen ${pathogen_genome} \
--reads "/vol/projects/bmikagos/Dual_RNA_seq_data/hela_Salmonella/fastq/*.fq.gz" \
--single_end --outdir ${out_dir_dualrnaseq} --run_salmon_alignment_based_mode --libtype ${libtype} \
--mapping_statistics --RNA_classes_to_replace_host ${RNA_classes_to_replace_host} \
--gene_feature_gff_to_create_transcriptome_pathogen "[gene, sRNA, tRNA, rRNA]" \
--gene_feature_gff_to_quantify_pathogen "[gene, sRNA, tRNA, rRNA]" \
--run_salmon_selective_alignment --run_star --run_htseq_uniquely_mapped
```


## 7. Read mapping and quantification

The nf-core/dualrnaseq pipeline provides three strategies to map and quantify your dual RNA-seq data. 

1) The first approach involves **Salmon with Selective alignment** which performs mapping using Selective alignment algorith and quantifies the reads. ```--run_salmon_selective_alignment```

2) The second approach involves **Salmon with alignment-based mode** which utalises aligned reads from STAR to quantify input reads. ```--run_salmon_alignment_based_mode```

3) The Third strategy utilizes alignment-based mapping executed with **STAR**, counting unique aligned reads with **HTSeq**. ```--run_star``` and ```--run_htseq_uniquely_mapped```


### 1) Salmon - Selective-Alignment

Salmon is a transcriptome-based mapping tool that performes both mapping and quantification. In the first phase it performes indexing of reference transcripts (pathogen transcripts are defined as gene or CDS features). A chimeric transcriptome of host and pathogen files is then created, utalising the associated fasta and gff files. Then coordinates of gene features are extracted from the host and pathogen annotation files.     

Salmon is a transcriptome-based mapping method. To avoid spurious mapping of reads that originate from unannotated locus with sequences similar to annotatated transcripts, a decoy-aware transcriptome is created and incorporated into the index. In the pipeline the decoy sequence is created from both pathogen and host entire genomes (which are both required).

> Note: **Selective-Alignment** is an improvement to the original alignment-free approach (also called quasi-mapping). In Selective-alignment, the best transcript for a read from a set of mappings is selected based on the alignment-based score instead of the the longest exact match - which increases the accuracy of the tool. See [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html) 
 


### 2) Salmon - quantification in alignment-based mode

In this [mode](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-alignment-based-mode), Salmon performs quantification utilizing an aligned BAM file. In the nf-core/dualrnaseq pipeline, the alignment file is generated with STAR. The first step involves indexing of a chimeric genome, created from the host and pathogen genome fasta files. Next, STAR performs an alignment, but for the purpose of the Salmon, it generates alignments translated into transcript coordinates. To learn more on this behavior, please see ```Output in transcript coordinates``` in [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) 

> Note: There are numerous STAR-based flags that can be modified - please see the full list of parameters.

> Salmon performs quantification based on a reference transcriptome. It is recommended to allow the pipeline to create a transcriptome using the provided genome fasta files and annotation gff files.

> To create the chimeric gff file used by STAR, and the chimeric transcriptome used by Salmon, please make sure the following parameters are defined: ```--gene_feature_gff_to_create_transcriptome_host```,  ```--gene_feature_gff_to_create_transcriptome_pathogen```, ```--gene_feature_gff_to_create_transcriptome_pathogen), ```--gene_attribute_gff_to_create_transcriptome_host```, ```--gene_attribute_gff_to_create_transcriptome_pathogen```    

> To quantify the alignments, please specify the parameters ```--libtype``` and ```--incompatPrior```. 


### 3) STAR - alignment-based genome mapping

STAR is a splice-aware alignment tool which aligns reads to a reference genome. In the nf-core/dualrnaseq pipeline, STAR generates an index of a chimeric genome, then identifies and maps spliced alignments across splice junctions. Therefore, the paths to host and pathogen genomes and annotative files must be provided either through iGenomes or directly using the appropriate parameters:  ```--genome_host```, ```--genome_pathogen```, ```--gff_host``` and ```--gff_pathogen```.



## 8. Output files
<br><br>
Regan to complete
<br><br>


## 9. Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

### AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

```--awsqueue```  The JobQueue that you intend to use on AWS Batch.

```--awsregion``` The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

```--awscli```  The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.





## 10. Maping statistics

To summarise the mapping statistics including total mapped reads, unmapped reads, host-specifc and pathogen-specific mapped reads.

```--mapping_statistics``` Default = True  

calc:
Count_total_reads
Count_total_trimmed_reads
Scatterplots

salmon - extract_processed_reads , unmapped

plot_salmon_mapping_stats_host_pathogen,
RNA class statistics - host, RNA_classes_to_replace.csv , def "data/RNA_classes_to_replace.csv" host
combined, each sample
