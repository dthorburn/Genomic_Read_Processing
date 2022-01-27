#!/usr/bin/env nextflow

// Pipeline developed for trimming and mapping genomic reads using TrimGalore or Trimmomatic and BWA-MEM. 
// Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
// Date last modified: 27/01/2022

def helpMessage() {
  log.info """
        Usage:
          You'll first need to update the paths and config file to reflect your environment and ensure you are in the same directory as the scripts, then:
          qsub NF_Mapping.sh
          
          If you require more advanced trimming options, you can skip the trimming steps and place trimmed gzipped fastqs into the 03_Trimmed directory and run:
          nextflow run Mapping.nf -c nextflow.config --profile imperial --Skip Trim 1

          The pipeline expects paired-end gzipped fastqc files that can be detected with the regex "*_R{1,2}*q.gz". 
          To check, use "ls -1 /path/to/02_Raw_Reads/*_R{1,2}*q.gz"
          
          Directory Structure:
            /Project_dir/                                                Project Directory - Exectute scripts from here
              | - Nextflow_Submit.sh                                     Pipeline submission script
              | - Mapping.nf                                             Nextflow script
              | - nextflow.config                                        Nextflow config - Update to reflect environment and computational requirements
              | - 01_FastQC/                                             
              | - 02_Raw_Reads/                                          Place all raw paired end read in this directory
                    | - 01_FastQC/                                       Optional post-trimming FastQC directory
              | - 03_Trimmed/
              | - 04_Mapped/
          
          Optional arguments:
            --help                                                       Show this message
            --init                                                       To be run first and only once - sets up conda environments
            --Skip_Trim                                                  Skips trimming step
            --FastQC                                                     Runs FastQC after trimming alongside mapping (Cannot be used with --Skip_Trim; off by default)
            --Skip_IndexRef                                              Skips the index reference step. Reference genome and bwa index files need to be in the same directory.    
            --mode trim_galore                                           Choice of which trimming software (trim_galore/trimmomatic; default: trim_galore)
           """
}

def versionMessage() {
  log.info  """
            Versions:
              pip/21.2.2
              conda/4.10.3
              samtools/1.2
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

// Set up conda environments if modules not working
if (params.init) {
  echo true

  params.Skip_Map = true
  params.Skip_Trim = true
  params.Skip_IndexRef = true
  
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

params.publishDir = './'
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

    beforeScript 'module load anaconda3/personal; source activate TrimGalore; module load samtools/1.2'
    
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
  Channel
    .fromFilePairs("${params.publishDir}/02_Raw_Reads/*_R{1,2}*q.gz")
    .ifEmpty { error "Cannot find any fastq files in ${params.publishDir}/02_Raw_Reads" }
    .set { raw_fastqs }
    
  process Trimming {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

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
    tuple val(sampleID), path("*R{1,2}*.fq.gz") into trimmed_fastqs, for_qc
    
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
          ${sampleID}_trimmed_R1_001.fq.gz ./01_Orphaned/${sampleID}_orphaned_R1_001.fastq.gz \\
          ${sampleID}_trimmed_R2_001.fq.gz ./01_Orphaned/${sampleID}_orphaned_R2_001.fastq.gz \\
          HEADCROP:${params.TG_Clip_R1} TRAILING:${params.TG_qual}
        """
    else
        error "Invalid alignment mode: ${params.mode}, use either trim_galore or trimmomatic"
  } 
}

// FastQC post-trimming
if( params.FastQC ) {
  process FastQC {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=8:mem=12gb -lwalltime=2:00:00"

    publishDir(
      path: "${params.publishDir}/02_Trimmed/01_FastQC",
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
      .fromFilePairs("${params.publishDir}/03_Trimmed/*_R{1,2}*q.gz")
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
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

    publishDir(
      path: "${params.publishDir}/04_Mapped",
      mode: 'copy',
    )
    
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BWA_threads}:mem=${params.BWA_memory} -lwalltime=${params.BWA_walltime}:00:00"
    
    input:
    path ref_genome
    path ref_index from ref_ch
    set trimmedID, path(trimmed_reads) from trimmed_fastqs
    
    output:
    // tuple val( "$trimmedID" ), path("${trimmedID}.bam"), path("${trimmedID}.bam.bai") into bams_ch // if you were to continue to develop this pipeline. 
    path("${trimmedID}.bam")
    path("${trimmedID}.bam.bai")
    
    beforeScript 'module load anaconda3/personal; source activate TrimGalore; module load samtools/1.2'
    
    script:
    def (read1, read2) = trimmed_reads
    // You have to use a single \ to escape $ defined in the bash script otherwise nextflow thinks they are nf variables and will throw an error. 
    """
    header=`zcat ${trimmedID}_R1_001.fastq.gz | head -n 1`
    RG_ID=`echo \$header | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/_/g'`
    RG_PU=`echo \$header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g'`

    ls

    bwa-mem2 mem -M -t ${params.BWA_threads} \\
      -R "@RG\\tID:\${RG_ID}\\tPU:\${RG_PU}\\tSM:${trimmedID}\\tLB:${params.BWA_RG_LB}\\tPL:${params.BWA_RG_PL}" \\
      ${ref_genome} ${read1} ${read2} > ${trimmedID}.sam
    samtools sort -o ./${trimmedID}.bam ${trimmedID}.sam
    samtools index ${trimmedID}.bam
    """
  }
}
