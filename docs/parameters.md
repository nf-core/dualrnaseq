# nf-core/dualrnaseq: Parameters

## Table of contents

1. [Nextflow](#1-nextflow)
2. [Genome references and annotation](#2-genome-references-and-annotation)
3. [Input sequence reads](#3-input-sequence-reads)
4. [FastQC and Adapter trimming](#4-fastqc-and-adapter-trimming)
5. [Salmon - general](#5-salmon---general)
6. [Salmon - selective alignment](#6-salmon---selective-alignment)
7. [STAR and Salmon - alignment based mode](#7-star-and-salmon---alignment-based-mode)
8. [HTSeq](#8-htseq)
9. [RNA mapping statistics](#9-rna-mapping-statistics)
10. [Reports](#10-reports)
11. [Other](#11-other)

### General comment

> All of the parameters listed here can be found in either the main configuration file `nextflow.config`, `base.config` or genome specific files such as `igenomes.conf` or  `genomes.conf`. Alternatively, each parameter can be specified by the user when they require adjustments to the default settings. The format for parameters is either a flag telling the pipeline to run something, such as `--run_STAR`, or to specify a particular value `--max_cpus 16`, string `--outWigStrand "Stranded"` or file `--outdir /path_to_file/file`.

## 1. Nextflow

### 1. Nextflow

> Note: most of the core Nextflow parameters only require a single hyphen

#### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
This is used in the MultiQC report and in the summary HTML / e-mail.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results where the inputs are the same, continuing from where it ran to previously.
Alternatively, you can supply a run name to resume a specific run: `-resume [run-name]`.
  
If unsure of the run name, you can use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file.

Note: you can use this to override pipeline defaults.

#### `--custom_config_version`

Provide git commit id for custom institutional configs hosted at `nf-core/configs`.

This was implemented for reproducibility purposes. Default: `master`

```bash
## Download and use config file with following git commit id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

#### `--custom_config_base`

If you're running offline, nextflow will be unable to fetch the institutional config files from the internet. If required, files should be downloaded when an active internet connection is available, then pointing Nextflow to those files with this option.

For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note: to make this process easier, the nf-core/tools helper package has a `download` command to download all required pipeline files + singularity containers + institutional configs in one go.

#### `--max_memory 128.GB`

Use to set the maximum memory for each process.

#### `--max_time 240.h`

Use to set the max run time for each process.

#### `--max_cpus 16`

Use to set the max number of CPUs for each process.

#### `--outdir`

Where results will be saved (should be enclosed by quotation marks `"..."`) .

### 2. Genome references and annotation

By default, Nextflow will locate and use the file `genomes.config` or `igenomes.config` to associate predefined genome references and annotations. If wishing to specify these manually, the following commands can be used.
> Note: we always recommend adding associated genome-based files to configuration files - to avoid clashes between user defined parameters and those supplied in configuration files.

#### `--igenomes_ignore`

To ignore all genome-based references in the igenomes configuration file

The following nine parameters are all set as `False`. If specified, the folder/file should be enclosed by quotations `"..."`.

#### `--fasta_host` 'path to file'

#### `--fasta_pathogen`

#### `--gff_host`

#### `--gff_host_tRNA`

#### `--gff_pathogen`

#### `--read_transcriptome_fasta_host_from_file`

#### `--read_transcriptome_fasta_pathogen_from_file`

#### `--transcriptome_host`

#### `--transcriptome_pathogen`

### 3. Input sequence reads

#### `--reads`

Input files can be read as either .fastq or .fastq.gz. They should be named descriptively without spaces and special characters (such as : and @), with the corresponding replicate (if any) appended at the end. The best practise for this pipeline is to use underscores to separate different experimental conditions.

Please note the following requirements:

* The path must be enclosed in quotes
* The path must have at least one `*` wildcard character
* When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.
* If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

> Note: by default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when launched. For example: `--single_end --reads '*.fastq'`

It is not possible to run a mixture of single-end and paired-end files in one run.

### 4. FastQC and adapter trimming

#### `--skipFastqc False`

An option to not run FastQC

> Note: Perhaps using BBduck would be easier - as it has an adaptor file built in with common methods including TruSeq etc

To remove adapter sequences that were introduced during the library preparation the pipeline utilizes cutadapt.
To learn more on cutadapt and its parameters visit the [`cutadapt documentation.`](https://cutadapt.readthedocs.io/en/stable/guide.html)

By default, the pipeline trims Illumina TruSeq adapters. See [`Illumina TruSeq.`](https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq)

#### `--skipTrimming False`

Will skip the trimming stage

#### `--a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"`

For single-end reads as well as the first reads of paired-end data, adapter sequence can be specified with `--a` flag. For more information, see [`adapter-types.`](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)

#### `--A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"`

For paired-end data, the adapter sequence for the second reads can be defined here. For more information, see [`trimming paired-end reads.`](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads)

#### `--quality_cutoff 10` or `--quality-cutoff 10,15`

Cutadapt can also remove low-quality read ends. By default, the 3’ end of each read is trimmed using a cutoff of 10. For more information on cutoff values, see [`quality trimming.`](https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming)

If you specify two comma-separated cutoffs, the first value represents the 5’ cutoff, and the second one the 3’ cutoff.

### 5. Salmon - general

#### `--libtype SF`

To define the type of sequencing library of your data.

To learn more on library types available in Salmon, please read [_`What’s this LIBTYPE?`_](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype)

#### `--kmer_length 21`

To define the k-mer length (`-k` parameter in Salmon, see [`preparing transcriptome indices`](https://salmon.readthedocs.io/en/latest/salmon.html?highlight=index#preparing-transcriptome-indices-mapping-based-mode)).

#### `--writeUnmappedNames True`

By default the pipeline saves names of unmapped reads. You can learn more about this option in [`Salmon documentation`](https://salmon.readthedocs.io/en/latest/salmon.html#writeunmappednames). If you don't want to keep this option, set the `--writeUnmappedNames` flag to false.

#### `--softclipOverhangs True`

By default, the pipeline allows soft-clipping of reads.

_"Soft-clipping allows reads that overhang the beginning or ends of the transcript. In this case, the overhanging section of the read will simply be unaligned, and will not contribute or detract from the alignment score"_.
If it is set to `False`, the end-to-end alignment of the entire read is forced, so that the occurance of any overhangings may affect the alignment score.

#### `--incompatPrior 0.0`

By default, this is set to `0.0`, to ensure that only mappings or alignments that are compatible with the specified library type are considered by Salmon. You can find more information on this parameter in the [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#incompatprior)

#### `--dumpEq True`

By default, to save the equivalence classes and their counts this option is set to `True`. See [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#dumpeq) for more information.

#### `--writeMappings False`

If set to `True`, the pipeline will create a `mapping.sam` file containing mapping information. To learn more on this option, please view the [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#writemappings)

### 6. Salmon - Selective alignment

#### `--run_salmon_selective_alignment False`

To run Salmon with Selective alignment

#### `--gene_feature_gff_to_create_transcriptome_host ["exon", "tRNA"]`

The pipeline uses gene features from the 3rd column of the host annotative file (gff/gtf) to extract the coordinates of transcripts to be quantified.

By default, the pipeline uses `exon` from the `--gff_host` file and `tRNA` from the `--gff_host_tRNA` file.

#### `--gene_feature_gff_to_create_transcriptome_pathogen ["gene", "sRNA", "tRNA", "rRNA"]`

The pipeline uses gene features from the 3rd column of the pathogen annotative fikle (gff/gtf) to extract the coordinates of transcripts to be quantified.

By default, the pipeline uses features as `gene`, `sRNA`, `tRNA` and `rRNA` from the `--gff_pathogen` file.

#### `--gene_attribute_gff_to_create_transcriptome_host ["transcript_id"]`

This flag defines the gene attribute from the 9th column of the host annotative (gff/gtf) file, where the transcript names are extracted.

By default, the pipeline extracts `transcript_id` from the `--gff_host` file.

#### `--gene_attribute_gff_to_create_transcriptome_pathogen ["locus_tag"]`

This flag defines the gene attribute from the 9th column of the pathogen annotative (gff/gtf) file, where transcript, genes or CDS regions are extracted.

By default, the pipeline extracts `locus_tag` from the `--gff_pathogen` file.

### 7. STAR and Salmon - alignment based mode

#### `--run_salmon_alignment_based_mode FALSE`

Option to run Salmn in alignment mode

#### `--run_star False`

Option to run STAR

#### `--outWigType "None"`

Used to generate signal outputs, such as "wiggle" and "bedGraph". To view all available signal types, please see the [`STAR documentation`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

By default, the pipeline does not generate any of these files.

#### `--outWigStrand "Stranded"`

Options are `Stranded` or `Unstranded` when defining the strandedness of wiggle/bedGraph output. See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--outSAMunmapped "Within"`

By default, the pipeline saves unmapped reads within the main BAM file. If you want to switch off this option, set the `--outSAMunmapped` flag to `None`.  See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more details.

For paired-end reads, the `KeepPairs` parameter will record the unmapped mates for each alignment, and will keep it adjacent to its mapped read (only affects multi-mapping reads).

#### `--outSAMattributes "Standard"`

To specify the attributes of the output BAm file. The default value is `Standard`, but there are a range of options if needed. Please see the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for the full list.

By default, the pipeline uses `Standard` option to keep NH HI AS nM SAM attributes.

#### `--outFilterMultimapNmax 20`

To specify the maximum number of loci a read is allowed to map to.
By default, this  option is set to 20 in the pipeline (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--outFilterType "BySJout"`

By default, the pipeline keeps reads containing junctions that passed filtering into the file `SJ.out.tab`. This option reduces the number of ”spurious” junctions. (ENCODE standard options for long RNA-seq pipeline). You can read more about the flag and its options in the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

#### `--quantTranscriptomeBan "Singleend"`

The nf-core/dualrnaseq pipeline runs STAR to generate a transcriptomic alignments. By default, it allows for insertions, deletions and soft-clips (`Singleend` option). To prohibit this behaviour, please specify `IndelSoftclipSingleend`. See the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignSJoverhangMin 8`

The number of minimum overhang for unannotated junctions can be changed here. By default, the pipeline uses 8. (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignSJDBoverhangMin 1`

The number of minimum overhang for annotated junctions can be changed here. See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--outFilterMismatchNmax 999`

To define a threshold for the number of mismatches to be allowed. By default, the pipeline uses a large number `999` to switch this filter off. (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--outFilterMismatchNoverReadLmax 0.04`

Here, you can define a threshold for a ratio of mismatches to *read* length. The alignment will be considered if the ratio is less than or equal to this value. For 2x100b, max number of mismatches is 0.04x200=8 for paired-end reads. (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignIntronMin 20`

By default, the nf-core dualrnaseq pipeline uses 20 as a minimum intron length. If the genomic gap is smaller than this value, it is considered as a deletion.
(ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignIntronMax 1000000`

The maximum intron length is set to 1,000,000 (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignMatesGapMax 1000000`

The maximum genomic distance between mates is 1,000,000 (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

### 8. HTSeq

#### `--run_htseq_uniquely_mapped False`

Used to run HTSeq-count and extract uniquely mapped reads from both the host and pathogen

#### `--stranded "yes"`

A parameter for the library type. Options include `"yes"` or `"no"`.

#### Host

#### `--gene_feature_gff_to_quantify_host ["exon","tRNA"]`

#### `--host_gff_atribute "gene_id"`

#### Pathogen

#### `--gene_feature_gff_to_quantify_pathogen ["gene", "sRNA", "tRNA", "rRNA"]`

#### `--pathogen_gff_atribute "locus_tag"`

The four parameters above are used to extract gene features from both the host and pathogen. These values may need to be changed, especially for the pathogen, as many different names exist, such as `ID`, `Gene`, `Name`, `locus_tag` etc etc.

A good idea is to view the accompanying annotative file and examine the fields within.

> Note: If a `tRNA.gff` file is included, it is assumed that it has the same gene atribute as the annotative file (gff/gtf), i.e. `gene_id`

### 9. RNA mapping statistics

#### `--mapping_statistics False`

Option to generate mapping statistics. This will create the following:

* Count the total number of reads before and after trimming
* Scatterplots comparing all replicates (separate for both host and pathogen reads)
* Plots of the % of mapped/quantified reads
* Plots of RNA-class statistics (as many types can be identified, the parameter below `--RNA_classes_to_replace_host` can help to summarise these)

#### `--RNA_classes_to_replace_host "$baseDir/data/RNA_classes_to_replace.csv"`

Located within the `data/` folder of dualrnaseq, this tab delimited file contains headers which groups similar types of RNA classes together. This helps to keep the RNA-class names simplified for plotting purposes.

Initially, the user can run the pipeline without this table (or remove the 'others' column because they may be interested in scRNAs). Depending on the requirements, the user can decide which types should be included/excluded or grouped together.  

### 10. Reports

#### `--email`

Set this parameter with your e-mail address to get a summary e-mail with details of the run when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

#### `--email_on_fail`

This works the same as `--email`, except emails are only sent if the workflow is not successful.

#### `--max_multiqc_email_size 25.MB`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached.

#### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

#### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

#### `--multiqc_config`

Specify path to a custom MultiQC configuration file.

### 11. Other

#### AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

#### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

#### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

#### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.
