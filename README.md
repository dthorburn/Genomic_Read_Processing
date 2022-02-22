# Genomic Read Mapping and Variant Calling
## Overview
These pipelines were developed to firstly process genomic reads through the trimming and mapping stages, then call variants using [GATK](https://gatk.broadinstitute.org/hc/en-us). The pipelines were developed using [Nextfow](https://www.nextflow.io/) version 20.10.0 (version on Imperial HPC; date: 25/01/2022). 

## Step 1: Genomic Read Trimming and Mapping
You should assess the reads using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to get suitable trimming paramaters. The quality of the reads I've used so far in the Crisanti lab only need simple trimming (adapter and random primer removal, and a simple trailing quality trim), so if you need more complex paramaters please contact me and I can update the pipeline to better reflect these needs.

#### Trimming 
There are currently two options for trimming: [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

#### Mapping
This is conducted using [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2). If you require a different genomic mapping tool, please let me know and I can add more options to the tool for versitility. 

### Usage
To use this pipeline follow these instructions:

  1. Run FastQC to identify appropriate trimming paramaters.
  2. Clone this repository using `git clone` and move the `Mapping.*` scripts into the project directory or download into the project working directory. I would recommend using a scratch directory such as the `ephemeral` Imperial HPC directory due to the generation of numerous sizeable files.
  3. Place all raw paired-end reads into a directory called `02_Raw_Reads`, or place already trimmed reads into the `03_Trimmed` directory and use the `--Skip_Trim` option in the `Mapping.sh` file.
  4. Update the `Mapping.config` file appropriately. The paramaters that are required to be updated are highlighted.
  5. Create the conda environment:
```
module load anaconda3/personal
conda create -n TrimGalore
source activate TrimGalore
conda install -c bioconda trim-galore
conda install -c bioconda bwa-mem2
```
  6. Submit the pipeline using `qsub Mapping.sh`


Below is the help message from `Mapping.nf`:
```
Usage:
  You'll first need to update the paths and config file to reflect your environment and ensure you are in the same directory as the scripts, then:
  qsub NF_Mapping.sh
  
  If you require more advanced trimming options, you can skip the trimming steps and place trimmed gzipped fastqs into the 03_Trimmed directory and run:
  nextflow run Mapping.nf -c nextflow.config --profile imperial --Skip Trim 1

  The pipeline expects paired-end fastqc files that can be detected with the glob "*_{R1,R2,1,2}{.fastq.gz,.fq.gz,.fastq,.fq,_001.fastq.gz,_001.fq.gz,_001.fastq,_001.fq}".
  If you are supplying trimmed reads, please ensure the name fits this a pattern like this sampleID_R1.fastq.gz or sampleID_1.fastq.gz
  
  Directory Structure:
    /Project_dir/                                                 Project Directory - Exectute scripts from here
      | - Mapping.sh                                              Pipeline coordinator submission script
      | - Mapping.nf                                              Nextflow script
      | - Mapping.config                                          Nextflow config - Update to reflect environment and computational requirements
      | - 01_FastQC/                                             
      | - 02_Raw_Reads/                                           Place all raw paired end read in this directory
            | - 01_FastQC/                                        Optional post-trimming FastQC directory
      | - 03_Trimmed/
      | - 04_Mapped/
  
  Optional arguments:
    --help                                                        Show this message
    --version                                                     See versions used to develop pipeline
    --Skip_Trim                                                   Skips trimming step
    --FastQC                                                      Runs FastQC after trimming alongside mapping (Cannot be used with --Skip_Trim; off by default)
    --Skip_IndexRef                                               Skips the index reference step. Reference genome and bwa index files need to be in the same directory.    
    --mode trim_galore                                            Choice of which trimming software (trim_galore/trimmomatic; default: trim_galore)
```
#### Planned Updates

Permit users to provide additional trimming paramaters like Trim Galore's `fastqc_args` paramater.

## Step 2: Variant Calling
The `GATK_Variants_Call.nf` pipeline is set up to call variants in either genomic DNA reads or in RNA-seq datasets. The pipelines were established using [GATK best practises](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). Ideally, reads will be mapped using the `Mapping.nf` pipeline for genomic DNA reads. This pipeline was developed for use on Imperial College's HPC. 

### Mode RNAseq
![a](https://drive.google.com/uc?id=1yLptERtjtDzx36vqAPZCMJgvbhASdsjq)
### Mode DNAseq
![a](https://drive.google.com/uc?id=1HKtzOeobgOVjCXEUE0-5378ocBz6Age7)

### Prerequisites on Imperial HPC

1. Bams are sorted and indexed. 
2. Index reference genome with `CreateSequenceDictionary` in [picard](https://gatk.broadinstitute.org/hc/en-us/articles/4414602399643-CreateSequenceDictionary-Picard-) and `faidx` from [samtools](http://www.htslib.org/doc/samtools-faidx.html).
3. Input bam files are suitable named. The names of the samples will be taken from the bam file name, where ABC123.bam will be ABC123 in the vcf file.

### Hard Filtering

Hard Filtering is necessary when SNP panels are unavailable. But each dataset will need paramaters that best fit, hence the optional inclusion of the `DRNASnp_Parse.{R,sh}` scripts here. They parse unfiltered VCFs and simply need run the nextflow pipeline with the `--Skip_VF` option to generate raw SNP VCFs, then update the `Project_dir` variable in the `DRNASnp_Parse.sh` file and include the optional `conda install` argument below. 

### Usage
To use this pipeline follow these instructions:

1. Clone this repository into the project directory.
2. Create the conda environment
```
module load anaconda3/personal
conda create -n NF_GATK
source activate NF_GATK
conda install -c bioconda gatk4
conda install r-vcfR r-dplyr r-ggpubr r-ggplot2 r-stringr r-data.table ## Only necessary for optional hard filtering annotation
```
3. Add the required (and optional) arguments to the command. This pipeline expects sorted and indexed bam files. 
4. Add the correct path to the project directory.
5. Submit the Nextflow coordinator with `qsub GATK_Variant_Call.sh`. 

Please note that all bams within the input directory will be considered as belonging to the same dataset. Also, due to the long runtimes of some of the processes, the long node is requred to run the coordinator, which is limited to 1 per user. 

Below is the help message from `GATK_Variant_Call.nf`:
```
Usage:
  These pipelines were developed using the best practises workflows set out by GATK for RNA-seq reads and genomic reads 
  (germline variant discovery): https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows

  (DNAseq): BP -> HCG -> DBI -> GVCF -> SV -------> VF
  (RNAseq): BP -> HC  -> -------------> SV -> MV -> VF

  This pipeline expects sorted bam files and will treat all bams in the input directory as a single dataset. 

  To use, there are 3 steps:
  1. Update project directory path in GATK_Variant_Call.sh 
  2. Add required arguments listed below
  3. Submit pipeline coordinator using qsub GATK_Variant_Call.sh

  If you require available HPC jobs for alternative scripts lower job concurrency options. 

  Required arguments:
    --RefGen                                        Path to reference fasta. Usage '--RefGen /path/to/genome.fasta'
    --BamDir                                        Path to input bam directory. Required even if skipping BP step. 

  Optional arguments:
    -w                                              Path to nextflow working directory. (Default: ./work)
    --help                                          Show this message
    --version                                       Show versions used to develop pipeline
    --VC_mode                                       Variant calling modes (RNAseq or DNAseq; default is DNAseq). 
                                                    Usage '--VC_mode RNAseq'.
    --Chrom                                         User defined chromosome selection (Default: all). 
                                                    Usage '--Chrom "AgamP4_2R,AgamP4_3R"'. Selection must be comma 
                                                    delimited in quotes and match the names of the contigs in the 
                                                    fasta index file.
    --VF_Filts                                      User defined filters to pass SNPs (Default: QUAL < 30.0, MQ < 40.0, 
                                                    SOR > 3.0, QD < 2.0, FS > 60.0). Each filter must be followed by filter 
                                                    name (usage: '--VF_Filts "--filter-expression QUAL < 30.0" 
                                                    --filter-name FAIL_QUAL --filter-expression etc..."').
    --MD_args                                       Optional arguments for MarkDuplicates         (BP; Both)
    --HC_args                                       Optional arguments for HaplotypeCaller        (RNAseq)
    --HCG_args                                      Optional arguments for HaplotypeCaller-GVCF   (DNAseq)
    --DBI_args                                      Optional arguments for GenomicsDBImport       (DNAseq)
    --GVCF_args                                     Optional arguments for GenotypeGVCF           (DNAseq)
    --SV_args                                       Optional arguments for SelectVariants         (Both)
    --MV_args                                       Optional arguments for vcf-merge              (RNAseq)
    --VF_args                                       Optional arguments for VariantFiltration      (Both)

  Concurrency arguments:                            Imperial HPC only permits 50 jobs per user. These options limit the 
                                                    number of concurrent processes running per step. NB. Multiple 
                                                    processes can be running at the same time.
    --BP_Forks                                      Default: 24 (Both)
    --HC_Forks                                      Default: 24 (RNAseq)
    --MV_Forks                                      Default: 24 (RNAseq)
    --HCG_Forks                                     Default: 24 (DNAseq)
    --DBI_Forks                                     Default: 15 (DNAseq) - All HCG processes must complete before DBI begins.
    --GVCF_Forks                                    Default: 15 (DNAseq)
    --SV_Forks                                      Default: 10 (Both)
    --VF_Forks                                      Default: 10 (Both)

  Debugging arguments:
    -resume                                         Resumes pipeline once errors are resolved. Usage: '-resume curious_borg' 
                                                    when log file shows "Launching `GATK_Variant_Call.nf` [curious_borg]"
    --Skip_BP                                       Skips processing bams
    --Skip_HC                                       Skips RNAseq HaplotypeCaller
    --Skip_HCG                                      Skips DNAseq HaplotypeCaller
    --Skip_DBI                                      Skips DNAseq GenomicsDBImport
    --Skip_GVCF                                     Skips DNAseq GenotypeGVCFs
    --Skip_SV                                       Skips SelectVariants
    --Skip_MV                                       Skips Merging VCFs
    --Skip_VF                                       Skips VariantFiltration
    --ProcBamDir                                    Path to processed bam directory  - Use if providing processed files
    --HCDir                                         Path to processed individual vcf directory
    --MVDir                                         Path to merged vcf directory
    --DBIDir                                        Path to processed DBI directory
    --GVCFDir                                       Path to processed GVCF directory
    --SVDir                                         Path to merged selected vcf directory - emits SNPs and Indels separately
    --VFDir                                         Path to merged filtered vcf directory
    --BP_threads                                    Number of threads for each subprocess - swap BP for process any acronym to 
                                                    alter other processes. (i.e., DBI_walltime = 24) 
    --BP_memory                                     Number of Gb of memory for each subprocess 
    --BP_walltime                                   Number of hours for each subprocess (72 is maximum) 
```

