#!/usr/bin/env nextflow

/*
 * Pipeline developed for trimming and mapping genomic reads using TrimGalore or Trimmomatic and BWA-MEM2. 
 * Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
 * Date last modified: 24/02/2022
 */

def helpMessage() {
  log.info """
Usage:
  If you require more advanced trimming options, you can skip the trimming steps and provide trimmed gzipped fastqs using the --TrimDir command. 
  If you require available HPC jobs for alternative scripts lower job concurrency options. 

  Required arguments:
    --RefGen                                                      Path to reference genome. Usage: '--RefGen /path/to/genome.fasta'
    --InDir                                                       Path to directory with raw fastqs (Not required if providing trimmed fastq.gzs)
    --TrimDir                                                     Path to directory with trimmed fastqs (Not required if providing raw fastq.gzs)
  
  Optional arguments:
    --help                                                        Show this message
    --FastQC                                                      Runs FastQC after trimming alongside mapping (Cannot be used with --Skip_Trim; off by default)
    --mode trim_galore                                            Choice of which trimming software (trim_galore/trimmomatic; default: trim_galore)

  Concurrency arguments:
    --Trim_Forks                                                  Number of concurrent trimming jobs. Default: 20
    --Map_Forks                                                   Number of concurrent mapping jobs. Default: 24
    --QC_Forks                                                    Number of concurrent fastqc jobs. Default: 5

  Debugging arguments:
    --Skip_Trim                                                   Skips trimming step.
    --Skip_IndexRef                                               Skips the index reference step. Reference genome and bwa index files need to be in the same directory.    
    --Skip_Map                                                    Skips the mapping step.
  """
}

def versionMessage() {
  log.info  """
            For repeatability, here are the versions I used to construct the pipeline:
              pip/21.2.2
              conda/4.10.3
              samtools/1.2
              fastqc/0.11.9
              cutadapt/1.18
              bwa-mem2/2.2.1
              trimmomatic/0.36
            """  
}

log.info """
==============================================================================================================================
                                        Mapping Reads Pipeline v2
==============================================================================================================================

Reference     : ${params.RefGen}
Trim Mode     : ${params.mode}
Input         : ${params.InDir}
Mapped        : ${PWD}/02_Mapped/

==============================================================================================================================
"""

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
if(!params.RefGen) {
  log.info"""
ERROR: No reference genome path provided! --RefGen /path/to/genome.fasta
==============================================================================================================================
  """
  helpMessage()
  exit 0
}
if(params.TrimDir == "./01_Trimmed"){
  if(!params.InDir) {
    log.info"""
ERROR: No input tarball provided! --InDir /path/to/raw_fastqs.gz
==============================================================================================================================
    """
    helpMessage()
    exit 0
  }
} else {
  if(!params.TrimDir) {
    log.info"""
ERROR: No trimmed fastq directory path provided! --TrimDir /path/to/trimmed_fastqs.gz
==============================================================================================================================
    """
    helpMessage()
    exit 0
  }
}


                                                            // =========================================================
                                                            // Setting the value channels (can be read unlimited times)
                                                            // =========================================================

ref_genome = file( params.refGen, checkIfExists: true )
ref_dir = ref_genome.getParent()
                                                            // =========================================================
                                                            // Step 1: Indexing Reference Genome
                                                            // =========================================================
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
    ## This is just to stop the pipeline. 
    exit 1
    """
  }
}

if( params.Skip_Trim == false ) {
  // use .flatten TRUE to permit emission as a single block rather than individial files.
  // File pairs are emitted as a tuple with this kind of structure [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
  // This glob pattern is terrible, and could easily lead to mistakes, but nextflow doesn't use the ?() syntax
  Channel
    .fromFilePairs("${params.InDir}/*_{R1,R2,1,2}{.fastq.gz,.fq.gz,.fastq,.fq,_001.fastq.gz,_001.fq.gz,_001.fastq,_001.fq}")
    .ifEmpty { error "Cannot find any fastq files in ${params.InDir}" }
    .set { raw_fastqs }
    
  process Trimming {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3
    maxForks params.Trim_Forks

    publishDir(
      path: "${params.TrimDir}",
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
          ${sampleID}_val_1.fq.gz ./01_Orphaned/${sampleID}_orphaned_R1.fastq.gz \\
          ${sampleID}_val_2.fq.gz ./01_Orphaned/${sampleID}_orphaned_R2.fastq.gz \\
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
      path: "${params.TrimDir}/01_FastQC",
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
      .fromFilePairs("${params.TrimDir}/*_{R1,R2,1,2}{.fastq.gz,.fq.gz,.fastq,.fq,_001.fastq.gz,_001.fq.gz,_001.fastq,_001.fq}")
      .ifEmpty { error "Cannot find any fastq files in ${params.TrimDir}" }
      .set { trimmed_fastqs }
  }
  // In case the reference index step is skipped, this should load the files into the working directory
  if( params.Skip_IndexRef ) {
    Channel
      .fromPath("${ref_dir}/*.f{asta,a,na}.*")
      .ifEmpty { error "Cannot find any reference index files in ${ref_dir}/" }
      .set { ref_ch }
  }

  process Mapping {
    // Sleeps 1 hour * attempt retries
    errorStrategy { sleep(Math.pow(2, task.attempt) * 3600 as long); return 'retry' }
    maxRetries 3
    maxForks params.Map_Forks

    publishDir(
      path: "${params.MapDir}",
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
    
    beforeScript 'module load anaconda3/personal; source activate TrimGalore; module load samtools/1.2'
    
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
