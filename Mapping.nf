#!/usr/bin/env nextflow

// Pipeline developed for trimming and mapping genomic reads using TrimGalore or Trimmomatic and BWA-MEM2. 
// Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
// Date last modified: 31/01/2022

def helpMessage() {
  log.info """
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
           """
}

def versionMessage() {
  log.info  """
            For repeatability, here are the versions I used to construct the pipeline:
              pip/21.2.2
              conda/4.10.3
              samtools/1.3.1
              fastqc/0.11.9
              cutadapt/1.18
              bwa-mem2/2.2.1
              trimmomatic/0.36
            """  
}

println "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nGenomic Read Mapping v0.1\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
// println "${PWD}, ${HOME}, ${PATH}"

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Show version information
if (params.version) {
    versionMessage()
    exit 0
}

// Set up conda environments if modules not working - deprecated
if (params.init) {
  echo true
  
  process Init {
    executor = 'local'
    beforeScript 'module load anaconda3/personal'
    
    script:
    """
    conda create -n TrimGalore
    source activate TrimGalore
    conda install -c bioconda trim-galore
    conda install -c bioconda bwa-mem2
    """
  }
}

params.publishDir = '.'
ref_genome = file( params.refGen )
ref_dir = ref_genome.getParent()

if( params.Skip_IndexRef == false ) {
  process IndexRef {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

    publishDir(
      path: "${ref_dir}/",
      mode: 'copy',
    )

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BWA_threads}:mem=${params.BWA_memory} -lwalltime=${params.BWA_walltime}:00:00"

    input:
    path ref_genome

    output:
    path "*fasta*" into ref_ch

    beforeScript 'module load anaconda3/personal; source activate TrimGalore; module load samtools/1.3.1'
    
    script:
    """
    bwa-mem2 index ${ref_genome}
    samtools faidx ${ref_genome}
    """
  }
}

if( params.Skip_Trim == false ) {
  // use .flatten TRUE to permit emission as a single block rather than individial files.
  // File pairs are emitted as a tuple with this kind of structure [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
  // This glob pattern is terrible, and could easily lead to mistakes, but nextflow doesn't use the ?() syntax
  Channel
    .fromFilePairs("${params.publishDir}/02_Raw_Reads/*_{R1,R2,1,2}{.fastq.gz,.fq.gz,.fastq,.fq,_001.fastq.gz,_001.fq.gz,_001.fastq,_001.fq}")
    .ifEmpty { error "Cannot find any fastq files in ${params.publishDir}/02_Raw_Reads" }
    .set { raw_fastqs }
    
  process Trimming {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3
    maxForks params.Trim_Forks

    publishDir(
      path: "${params.publishDir}/03_Trimmed",
      mode: 'copy',
    )
    
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.TG_threads}:mem=${params.TG_memory} -lwalltime=${params.TG_walltime}:00:00"
    
    input:
    // Because it's a tuple channel, the first entry, the common characters among samples will be the sampleID, which is the easiest to work with, but you still need to declare the other entries (I think). 
    set sampleID, path(reads) from raw_fastqs
    
    output:
    // outputs a structure like this: [sampleID, [Read1, Read2]]
    tuple val(sampleID), path("*_val_{1,2}.fq.gz") into trimmed_fastqs, for_qc
    
    // Unsure if you can do contidional beforeScript arguments so put environment loading into the script    
    beforeScript 'module load anaconda3/personal; source activate TrimGalore; module load trimmomatic/0.36'
      
    script:
    def (read1, read2) = reads
    
    if( params.mode == 'trim_galore' )
        """
        trim_galore --paired --gzip -j 4 \\
          --clip_R1 ${params.TG_Clip_R1} --clip_R2 ${params.TG_Clip_R2} \\
          --quality ${params.TG_qual} ${read1} ${read2}
        """
    else if( params.mode == 'trimmomatic' )
        """
        mdkir 01_Orphaned
        trimmomatic PE -validatePairs -threads ${params.TG_threads} ${read1} ${read2} \\
          ${sampleID}_trimmed_R1.fq.gz ./01_Orphaned/${sampleID}_orphaned_R1.fastq.gz \\
          ${sampleID}_trimmed_R2.fq.gz ./01_Orphaned/${sampleID}_orphaned_R2.fastq.gz \\
          HEADCROP:${params.TG_Clip_R1} TRAILING:${params.TG_qual}
        """
    else
        error "Invalid alignment mode: ${params.mode}, use either trim_galore or trimmomatic"
  } 
}

// FastQC post-trimming
if( params.FastQC ) {
  process FastQC {
    // Sleeps 1 hour * attempt retries
    errorStrategy { sleep(Math.pow(2, task.attempt) * 3600 as long); return 'retry' }
    maxRetries 3
    maxForks params.QC_Forks

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=8:mem=12gb -lwalltime=2:00:00"

    publishDir(
      path: "${params.publishDir}/03_Trimmed/01_FastQC",
      mode: 'move',
    )

    input:
    set sampleID, path(reads) from for_qc

    output:
    path "*.html"

    beforeScript 'module load fastqc/0.11.9'

    script:
    def (read1, read2) = reads
    """
    mkdir tmp
    fastqc -t 8 -d ./tmp/ ${read1}
    fastqc -t 8 -d ./tmp/ ${read2}
    """
  }
}

if( params.Skip_Map == false ) {
  
  // This is just in case trimming is being skipped, the input files will still be loaded in to a channel.
  if( params.Skip_Trim ) {
    Channel
      .fromFilePairs("${params.publishDir}/03_Trimmed/*_{R1,R2,1,2}{.fastq.gz,.fq.gz,.fastq,.fq,_001.fastq.gz,_001.fq.gz,_001.fastq,_001.fq}")
      .ifEmpty { error "Cannot find any fastq files in ${params.publishDir}/03_Trimmed" }
      .set { trimmed_fastqs }
  }
  // In case the reference index step is skipped, this should load the files into the working directory
  if( params.Skip_IndexRef ) {
    Channel
      .fromPath("${ref_dir}/*.fasta.*")
      .ifEmpty { error "Cannot find any reference index files in ${ref_dir}/" }
      .set { ref_ch }
  }

  process Mapping {
    // Sleeps 1 hour * attempt retries
    errorStrategy { sleep(Math.pow(2, task.attempt) * 3600 as long); return 'retry' }
    maxRetries 3
    maxForks params.Map_Forks

    publishDir(
      path: "${params.publishDir}/04_Mapped",
      mode: 'move',
    )
    
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BWA_threads}:mem=${params.BWA_memory} -lwalltime=${params.BWA_walltime}:00:00"
    
    input:
    path ref_genome
    path ref_index from ref_ch.collect()
    set trimmedID, path(trimmed_reads) from trimmed_fastqs
    
    output:
    // tuple val( "$trimmedID" ), path("${trimmedID}.bam"), path("${trimmedID}.bam.bai") into bams_ch // if you were to continue to develop this pipeline. 
    path("${trimmedID}.bam")
    path("${trimmedID}.bam.bai")
    
    beforeScript 'module load anaconda3/personal; source activate TrimGalore; module load samtools/1.3.1'
    
    script:
    def (read1, read2) = trimmed_reads
    // You have to use a single \ to escape $ defined in the bash script otherwise nextflow thinks they are nf variables and will throw an error. 
    """
    header=`zcat ${read1} | head -n 1`
    RG_ID=`echo \$header | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/_/g'`
    RG_PU=`echo \$header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g'`

    bwa-mem2 mem -M -t ${params.BWA_threads} \\
      -R "@RG\\tID:\${RG_ID}\\tPU:\${RG_PU}\\tSM:${trimmedID}\\tLB:${params.BWA_RG_LB}\\tPL:${params.BWA_RG_PL}" \\
      ${ref_genome} ${read1} ${read2} > ${trimmedID}.sam
    samtools sort -o ./${trimmedID}.bam ${trimmedID}.sam
    samtools index ${trimmedID}.bam
    """
  }
}
