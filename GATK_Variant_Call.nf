#!/usr/bin/env nextflow

/*  
 * 
 * 'GATK_Variant_Call.nf' is a variant calling pipeline paired with the 'Mapping.nf' pipeline.
 * 
 * Date last modified: 19/08/2022
 * 
 * Authors:
 * Miles Thorburn
 * Nace Kranjc
 * 
 */

def helpMessage() {
  log.info """
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
        --Chroms                                        User defined chromosome selection (Default: all). 
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
        --MV_args                                       Optional arguments for MergeVcfs              (RNAseq)
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

    ==============================================================================================================================
  """
}

def versionMessage() {
  log.info """
    For repeatability, here are the versions I used to construct the pipeline:
    gatk/4.2.4.1
    conda/4.10.3
    samtools/1.2
    vcftools/1.3.1
    bash/4.4.20
           """
}

log.info """
==============================================================================================================================
                                      GATK Variant Calling Nextflow Pipeline v1
==============================================================================================================================

Mode          : ${params.VC_mode}
Reference     : ${params.RefGen}
InBams        : ${params.BamDir}
RawVars       : ${PWD}/06_SelectVariants/
FiltSnps      : ${PWD}/07_Filtered_SNPs/

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

if(!params.BamDir) {
  log.info"""
ERROR: No bam path provided! --BamDir /path/to/bams/
==============================================================================================================================
  """
  helpMessage()
  exit 0
}

// ============================================
// Setting the value channels (can be read unlimited times)
// ============================================
ref_genome = file( params.RefGen, checkIfExists: true )
ref_dir    = ref_genome.getParent()
ref_name   = ref_genome.getBaseName()
ref_dict   = file( "${ref_dir}/${ref_name}.dict", checkIfExists: true )
ref_index  = file( "${ref_dir}/${ref_name}*.fai", checkIfExists: true ).first()

// ============================================
// Processing Bams - 1 step
// ============================================
if( params.Skip_BP == false ){
  Channel
    .fromPath("${params.BamDir}/*.bam")
    .ifEmpty { error "No bams found in ${params.BamDir}" }
    .map { file -> tuple(file.baseName, file) }
    .set { raw_bam_ch }

  Channel
    .fromPath("${params.BamDir}/*.bai")
    .ifEmpty { error "No bais found in ${params.BamDir}" }
    .map { file -> tuple(file.baseName.replaceAll(".bam", ""), file) } // Handles the cases where files are names SampleID.bam.bai
    .set { raw_bai_ch }

  raw_bam_ch
    .combine( raw_bai_ch, by: 0 )
    .set { all_bams }

  process GATK_BP {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.BP_Forks

    tag { SampleID }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BP_threads}:mem=${params.BP_memory}gb -lwalltime=${params.BP_walltime}:00:00"

    publishDir(
      path: "${params.ProcBamDir}",
      mode: 'copy',
    )

    input:
    tuple SampleID, path(bam), path(bai) from all_bams

    output:
    tuple val(SampleID), path("${SampleID}_BP.bam"), path("${SampleID}_BP.bam.bai") into processed_bams
    path("*.txt")

    beforeScript 'module load anaconda3/personal; source activate NF_GATK; module load samtools/1.2'

    script:
    if( params.VC_mode == "DNAseq" )
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.BP_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk MarkDuplicates \\
        --TMP_DIR tmp/ \\
        -I ${bam} \\
        -O ${SampleID}_BP.bam \\
        -M ${SampleID}_MD.txt ${params.MD_args}

      samtools index ${SampleID}_BP.bam
      """
    else if( params.VC_mode == "RNAseq")
      // May need to break this up
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.BP_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi

      taskset -c 0-\${n_slots} gatk MarkDuplicates \\
        --TMP_DIR tmp/ \\
        -I ${bam} \\
        -O ${SampleID}_MD.bam \\
        -M ${SampleID}_MD.txt ${params.MD_args}

      echo "Finished with MD for ${SampleID}:`date`"

      taskset -c 0-\${n_slots} gatk SplitNCigarReads \\
        --tmp-dir tmp/ \\
        -R ${ref_genome} \\
        -I ${SampleID}_MD.bam \\
        -O ${SampleID}_BP.bam

      samtools index ${SampleID}_BP.bam
      """ 
      else
        error "Invalid variant calling mode (${params.VC_mode}), use either DNAseq or RNAseq"
  }
}
                                                            // ============================================
                                                            // RNAseq Pipeline Start - 2 Steps
                                                            // ============================================
if( params.VC_mode == "RNAseq" ){
                                                            // ============================================
                                                            // RNAseq: Step 1 - HaplotypeCaller
                                                            // ============================================
  if( params.Skip_HC == false ){
    // Reading in data if first step is skipped
    if( params.Skip_BP ){
      Channel
        .fromPath("${params.BamDir}/*.bam")
        .ifEmpty { error "No bams found in ${params.ProcBamDir}" }
        .map { file -> tuple(file.baseName, file) }
        .set { proc_bam_ch }

      Channel
        .fromPath("${params.BamDir}/*.bai")
        .ifEmpty { error "No bais found in ${params.ProcBamDir}" }
        .map { file -> tuple(file.baseName.replaceAll(".bam", ""), file) }
        .set { proc_bai_ch }

      proc_bam_ch
        .combine( proc_bai_ch, by: 0 )
        .set { processed_bams }
    }
    // Setting up the chromosome channel
    if( params.Chroms == "" ){
      // Defaulting to using all chromosomes
      chromosomes_ch = Channel
                          .from("AgamP4_2L", "AgamP4_2R", "AgamP4_3L", "AgamP4_3R", "AgamP4_X", "AgamP4_Y_unplaced", "AgamP4_UNKN")
      println "No chromosomes specified, using all major chromosomes: AgamP4_2L, AgamP4_2R, AgamP4_3L, AgamP4_3R, AgamP4_X, AgamP4_Y_unplaced, AgamP4_UNKN"
    } else {
      // User option to choose which chromosome will be used
      // This worked with the following syntax nextflow run testing.nf --profile imperial --Chroms "AgamP4_3R,AgamP4_2L"
      chrs = params.Chroms.split(",")
      chromosomes_ch = Channel
                          .from(chrs)
      println "User defined chromosomes set: ${params.Chroms}"
    }

    process RNA_HC {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
      maxRetries 3
      maxForks params.HC_Forks

      tag { SampleID+"-"+chrom }

      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.HC_threads}:mem=${params.HC_memory}gb:mpiprocs=1:ompthreads=${params.HC_threads} -lwalltime=${params.HC_walltime}:00:00"

      publishDir(
        path: "${params.HCDir}",
        mode: 'copy',
      )

      input:
      each chrom from chromosomes_ch
      set SampleID, path(bam), path(bai) from processed_bams
      path ref_genome
      path ref_dict
      path ref_index

      output:
      //tuple chrom, path("${SampleID}-${chrom}.vcf") into MV_ch
      //path("${SampleID}-${chrom}.vcf.idx") into MV_idxs
      tuple sampleID, chrom, path("${SampleID}-${chrom}.vcf"), path("${SampleID}-${chrom}.vcf.idx") into SV_ch
      
      beforeScript 'module load anaconda3/personal; source activate NF_GATK'

      script:
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.HC_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk --java-options \"-Xmx${params.HCG_memory}G -XX:+UseParallelGC -XX:ParallelGCThreads=\${n_slots}\" HaplotypeCaller \\
        --tmp-dir tmp/ \\
        --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP \\
        --native-pair-hmm-threads \${n_slots} \\
        -R ${ref_genome} \\
        -I ${bam} \\
        -O ${SampleID}-${chrom}.vcf \\
        -L ${chrom} \\
        --standard-min-confidence-threshold-for-calling 20 ${params.HC_args}
      """
    }
  // Deprecated due to issues with merging variants
  //  MV_ch
  //    .groupTuple(by: 0)
  //    .set{ MV_chr_ch }
  }

// ============================================
// RNAseq: Step 2 - SelectVariants
// ============================================
  if( params.Skip_SV == false ){
    if( params.Skip_HC ){
      Channel
        .fromFilePairs("${params.HCDir}/*{vcf,idx}") { file -> file.name.replaceAll(/.vcf|.idx$/,'') }
        .ifEmpty { error "No vcfs files found in ${params.GVCFDir}" }
        .map { ID, files -> tuple(ID.replaceAll("-.*",""), ID.replaceAll("^.*?-",""), files[0], files[1])}
        .set { SV_ch }
    }
      
    
    process RNA_SV {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
      maxRetries 3
      maxForks params.SV_Forks

      tag { SampleID+"-"+chrom }

      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.SV_threads}:mem=${params.SV_memory}gb -lwalltime=${params.SV_walltime}:00:00"

      publishDir(
        path: "${params.SVDir}",
        mode: 'copy',
      )

      input:
      //each chrom from chromosomes_ch
      //set chrom, path(vcf), path(idx) from sv_ch
      tuple SampleID, chrom, path(vcf), path(idx) from SV_ch
      path ref_genome
      path ref_dict
      path ref_index

      output:
      tuple chrom, path("${chrom}-${SampleID}_SNPs.raw.vcf"), path("${chrom}-${SampleID}_SNPs.raw.vcf.idx") into MV_ch
      tuple chrom, path("${chrom}-${SampleID}_INDELs.raw.vcf"), path("${chrom}-${SampleID}_INDELs.raw.vcf.idx") into raw_indel_ch
      
      beforeScript 'module load anaconda3/personal; source activate NF_GATK'

      script:
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.SV_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk SelectVariants \\
        --tmp-dir tmp/ \\
        -R ${ref_genome} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        --select-type-to-include MNP \\
        -O ${chrom}-${SampleID}_SNPs.raw.vcf ${params.SV_args}

      taskset -c 0-\${n_slots} gatk SelectVariants \\
        --tmp-dir tmp/ \\
        -R ${ref_genome} \\
        -V ${vcf} \\
        --select-type-to-include INDEL \\
        -O ${chrom}-${SampleID}_INDELs.raw.vcf ${params.SV_args}
      """
    }
    MV_ch
      .groupTuple(by: 0)
      .set{ MV_chr_ch }
  }

// ============================================
// RNAseq: Step 3 - MergeVcfs
// ============================================
  if( params.Skip_MV == false ){
    if( params.Skip_SV ){
      Channel
        .fromFilePairs("${params.SVDir}/*_SNPs.raw.{vcf,idx,vcf.idx}") { file -> file.name.replaceAll(/.vcf|.idx$/,'') }
        .ifEmpty { error "No vcfs files found in ${params.SVDir}" }
        .map { ID, files -> tuple(ID.replaceAll("-.*",""), files[0], files[1])}
        .groupTuple(by: 0)
        .set { MV_chr_ch }

    // Deprecated for new pipeline order. 
    //  Channel
    //    .fromPath("${params.SVDir}/*.vcf")
    //    .ifEmpty { error "No vcfs files found in ${params.SVDir}" }
    //    .map { file -> tuple(file.baseName.replaceAll(".vcf", "").replaceAll("^.*?-",""), file) }
    //    .groupTuple(by: 0)
    //    .set { MV_chr_ch }

    //  Channel
    //    .fromPath("${params.SVDir}/*.idx")
    //    .ifEmpty { error "No index files found in ${params.SVDir}" }
    //    .set { MV_idxs }
    }
    
    process RNA_MV {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
      maxRetries 3
      maxForks params.MV_Forks

      tag { chrom }

      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.MV_threads}:mem=${params.MV_memory}gb -lwalltime=${params.MV_walltime}:00:00"

      publishDir(
        path: "${params.MVDir}",
        mode: 'copy',
      )

      input:
      //each chrom from chromosomes_ch
      tuple chrom, path(vcf), path(idx) from MV_chr_ch
      //path(idx) from MV_idxs.collect()

      output:
      tuple chrom, path("${chrom}.vcf"), path("${chrom}.vcf.idx") into vf_ch
      
      beforeScript 'module load vcftools/0.1.13; module load htslib/1.3.2; module load anaconda3/personal; source activate NF_GATK'

      script:
      // Handling variable numbers of files being included
      //def vcf_params = vcf.collect{ "-I $it" }.join(' ')
      def vcf_compress = vcf.collect{ "bgzip $it; tabix -p vcf ${it}.gz" }.join('\n')
      def vcf_params = vcf.collect{ "${it}.gz" }.join(' ')
      """
      ${vcf_compress}

      vcf-merge -c snps -d ${params.MV_args} ${vcf_params} > ${chrom}.vcf

      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.MV_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk IndexFeatureFile --tmp-dir tmp/ -I ${chrom}.vcf
      """
    }
  }

// ============================================
// DNAseq Pipeline Start
// ============================================
} else if( params.VC_mode == "DNAseq" ){

// ============================================
// DNAseq: Step 1 - HaplotypeCaller GVCF mode
// ============================================
  if( params.Skip_HCG == false ){
    if( params.Skip_BP ){
      Channel
        .fromPath("${params.ProcBamDir}/*.bam")
        .ifEmpty { error "No bams found in ${params.ProcBamDir}" }
        .map { file -> tuple(file.baseName, file) }
        .set { proc_bam_ch }

      Channel
        .fromPath("${params.ProcBamDir}/*.bai")
        .ifEmpty { error "No bais found in ${params.ProcBamDir}" }
        .map { file -> tuple(file.baseName.replaceAll(".bam", ""), file) } // Handles the cases where files are names SampleID.bam.bai
        .set { proc_bai_ch }

      proc_bam_ch
        .combine( proc_bai_ch, by: 0 )
        .set { processed_bams }
    }
    // Setting up the chromosome channel
    if( params.Chroms == "" ){
      // Defaulting to using all chromosomes
      chromosomes_ch = Channel
                          .from("AgamP4_2L", "AgamP4_2R", "AgamP4_3L", "AgamP4_3R", "AgamP4_X", "AgamP4_Y_unplaced", "AgamP4_UNKN")
      println "No chromosomes specified, using all major chromosomes: AgamP4_2L, AgamP4_2R, AgamP4_3L, AgamP4_3R, AgamP4_X, AgamP4_Y_unplaced, AgamP4_UNKN"
    } else {
      // User option to choose which chromosome will be used
      // This worked with the following syntax nextflow run testing.nf --profile imperial --Chroms "AgamP4_3R,AgamP4_2L"
      chrs = params.Chroms.split(",")
      chromosomes_ch = Channel
                          .from(chrs)
      println "User defined chromosomes set: ${params.Chroms}"
    }

    process DNA_HCG {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
      maxRetries 3
      maxForks params.HCG_Forks

      tag { SampleID+"-"+chrom }

      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.HCG_threads}:mem=${params.HCG_memory}gb:mpiprocs=1:ompthreads=${params.HCG_threads} -lwalltime=${params.HCG_walltime}:00:00"

      publishDir(
        path: "${params.HCDir}",
        mode: 'copy',
      )

      input:
      each chrom from chromosomes_ch
      tuple SampleID, path(bam), path(bai) from processed_bams
      path ref_genome
      path ref_dict
      path ref_index

      output:
      tuple chrom, path("${SampleID}-${chrom}.vcf") into HCG_ch
      path("${SampleID}-${chrom}.vcf.idx") into idx_ch
      
      beforeScript 'module load anaconda3/personal; source activate NF_GATK'

      script:
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.HCG_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk --java-options \"-Xmx${params.HCG_memory}G -XX:+UseParallelGC -XX:ParallelGCThreads=\${n_slots}\" HaplotypeCaller \\
        --tmp-dir tmp/ \\
        --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP \\
        --native-pair-hmm-threads \${n_slots} \\
        -ERC GVCF \\
        -L ${chrom} \\
        -R ${ref_genome} \\
        -I ${bam} \\
        -O ${SampleID}-${chrom}.vcf ${params.HCG_args}
      """
    }
    // Collecting and grouping output of this process by chromosome. 
    // This operator collects output, alongside the .collect in idx_ch below, this entire process needs to complete
    // before DBI starts. This could be optimsed, but will only decrease runtime by a few hours at most
    HCG_ch
      .groupTuple(by: 0)
      .set{ HCG_chr_ch }
  }

// ============================================
// DNAseq: Step 2 - GenomicsDBImport
// ============================================
  if( params.Skip_DBI == false ){
    if( params.Skip_HCG ){
      Channel
        .fromPath("${params.HCDir}/*.vcf")
        .ifEmpty { error "No vcfs found in ${params.HCDir}" }
        .map { file -> tuple(file.baseName.replaceAll(".vcf", "").replaceAll("^.*?-",""), file) }
        .groupTuple(by: 0)
        .set { HCG_chr_ch }

      Channel
        .fromPath("${params.HCDir}/*.idx")
        .ifEmpty { error "No index files found in ${params.HCDir}" }
        .set { idx_ch }

      // Deprecated for chromosome parallisation
      //vcf_ch
      //  .combine( idx_ch, by: 0 )
      //  .set { HCG_ch }
    }
    // Deprecated for chromosome parallisation in HCG
    // Combining the input channel. Important which comes first in how the value is added to the tuple, here it is "ID, vcf, idx, chrom"
    // .combine is needed to iterate input channel, otherwise it will simply merge the tuples
    // Example .view() of Parallel_chrom_GVCF_ch is [AE12A_S24, ./01_Input_Bams/AE12A_S24.bam, ./01_Input_Bams/AE12A_S24.bam.bai, AgamP4_2L]
    //GVCF_ch
    //  .combine( chromosomes_ch )
    //  .set { Parallel_chrom_GVCF_ch }
    
    process DNA_DBI {
      errorStrategy { sleep(600); return 'retry' }
      maxRetries 3
      maxForks params.DBI_Forks

      tag { chrom }

      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.DBI_threads}:mem=${params.DBI_memory}gb -lwalltime=${params.DBI_walltime}:00:00"

      publishDir(
        path: "${params.DBIDir}",
        mode: 'copy',
      )

      input:
      //each chrom from chromosomes_ch - now split by chrom in the above process
      tuple chrom, path(vcf) from HCG_chr_ch
      // A bit lazy, but this offers all index files to every instance. They're symlinks, so it's not limiting in any way. 
      path(idx) from idx_ch.collect()
      path ref_genome
      path ref_dict
      path ref_index

      output:
      tuple chrom, path("DB.${chrom}", type: 'dir') into DBI_ch
      
      beforeScript 'module load anaconda3/personal; source activate NF_GATK'

      script:
      // Handling variable numbers of files being included
      def vcf_params = vcf.collect{ "-V $it" }.join(' ')
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      DB_name="DB.${chrom}"

      n_slots=`expr ${params.DBI_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk GenomicsDBImport \\
        --tmp-dir tmp/ \\
        ${vcf_params} \\
        --genomicsdb-workspace-path \${DB_name} \\
        -L ${chrom} ${params.DBI_args}
      """
    }
  }

// ============================================
// DNAseq: Step 3 - GenotypeGVCF
// ============================================
  if( params.Skip_GVCF == false ){
    if( params.Skip_DBI ){
      Channel
        .fromPath("${params.DBIDir}/*", type: 'dir')
        .map { it -> tuple(it.Extension, it) }
        .set { DBI_ch }
    }

    process DNA_GVCF {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
      maxRetries 3
      maxForks params.GVCF_Forks

      tag { chrom }

      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.GVCF_threads}:mem=${params.GVCF_memory}gb -lwalltime=${params.GVCF_walltime}:00:00"

      publishDir(
        path: "${params.GVCFDir}",
        mode: 'copy',
      )

      input:
      tuple chrom, path(database) from DBI_ch
      path ref_genome
      path ref_dict
      path ref_index

      output:
      tuple chrom, path("${chrom}.vcf"), path("${chrom}.vcf.idx") into sv_ch
      
      beforeScript 'module load anaconda3/personal; source activate NF_GATK'

      script:
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.GVCF_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk GenotypeGVCFs \\
        --tmp-dir tmp/ \\
        -R ${ref_genome} \\
        -V gendb://${database} \\
        -O ${chrom}.vcf \\
        -L ${chrom} ${params.GVCF_args}
      """
    }
  }

// ============================================
// DNAseq: Step 4 - SelectVariants
// ============================================
  if( params.Skip_SV == false ){
    if( params.Skip_GVCF ){
      // DNAseq input
      Channel
        .fromFilePairs("${params.GVCFDir}/*{vcf,idx}") { file -> file.name.replaceAll(/.vcf|.idx$/,'') }
        .ifEmpty { error "No vcfs files found in ${params.GVCFDir}" }
        .map { ID, files -> tuple(ID, files[0], files[1])}
        .set { sv_ch }
    } //else if( params.Skip_MV ){
      // RNAseq input
      //Channel
      //  .fromFilePairs("${params.MVDir}/*{vcf,idx}") { file -> file.name.replaceAll(/.vcf|.idx$/,'') }
      //  .ifEmpty { error "No vcfs files found in ${params.MVDir}" }
      //  .map { ID, files -> tuple(ID, files[0], files[1])}
      //  .set { sv_ch }
    //}
    
    process DNA_SV {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
      maxRetries 3
      maxForks params.SV_Forks

      tag { chrom }

      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.SV_threads}:mem=${params.SV_memory}gb -lwalltime=${params.SV_walltime}:00:00"

      publishDir(
        path: "${params.SVDir}",
        mode: 'copy',
      )

      input:
      //each chrom from chromosomes_ch
      set chrom, path(vcf), path(idx) from sv_ch
      path ref_genome
      path ref_dict
      path ref_index

      output:
      tuple chrom, path("${chrom}_SNPs.raw.vcf"), path("${chrom}_SNPs.raw.vcf.idx") into vf_ch
      tuple chrom, path("${chrom}_INDELs.raw.vcf"), path("${chrom}_INDELs.raw.vcf.idx") into raw_indel_ch
      
      beforeScript 'module load anaconda3/personal; source activate NF_GATK'

      script:
      """
      if [ ! -d tmp ]; then mkdir tmp; fi
      n_slots=`expr ${params.SV_threads} / 2 - 3`
      if [ \$n_slots -le 0 ]; then n_slots=1; fi
      taskset -c 0-\${n_slots} gatk SelectVariants \\
        --tmp-dir tmp/ \\
        -R ${ref_genome} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        --select-type-to-include MNP \\
        -O ${chrom}_SNPs.raw.vcf ${params.SV_args}

      taskset -c 0-\${n_slots} gatk SelectVariants \\
        --tmp-dir tmp/ \\
        -R ${ref_genome} \\
        -V ${vcf} \\
        --select-type-to-include INDEL \\
        -O ${chrom}_INDELs.raw.vcf ${params.SV_args}
      """
    }
  }
}

// ============================================
// VariantFiltration
// ============================================
if( params.Skip_VF == false ){
  if( params.Skip_SV ){
    if("${params.VC_mode}" == "DNAseq"){ 
      Channel
        .fromFilePairs("${params.SVDir}/*_SNPs.raw.{vcf,idx,vcf.idx}") { file -> file.name.replaceAll(/.vcf|.idx$/,'') }
        .ifEmpty { error "No vcfs files found in ${params.SVDir}" }
        .map { ID, files -> tuple(ID.replaceAll("_SNPs.raw", ""), files[0], files[1])}
        .set { vf_ch }
    } else if("${params.VC_mode}" == "RNAseq"){
      Channel
        .fromFilePairs("${params.MVDir}/*_SNPs.raw.{vcf,idx,vcf.idx}") { file -> file.name.replaceAll(/.vcf|.idx$/,'') }
        .ifEmpty { error "No vcfs files found in ${params.MVDir}" }
        .map { ID, files -> tuple(ID.replaceAll("_SNPs.raw", ""), files[0], files[1])}
        .set { vf_ch }
    }
  }

  // Setting up the filter channel
  if( params.VF_Filts == "" ){
    // Defaulting to using all chromosomes
    filter_ch = Channel.value(" --filter-expression \"QUAL < 30.0\" --filter-name FAIL_QUAL --filter-expression \"MQ < 40.0\" --filter-name FAIL_MQ --filter-expression \"SOR > 3.00\" --filter-name FAIL_SOR --filter-expression \"QD < 2.0\" --filter-name FAIL_QD --filter-expression \"FS > 60.0\" --filter-name FAIL_FS")
    println "No SNP filters specified, using defaults to fail variants: QUAL < 30.0, MQ < 40.0, SOR > 3.0, QD < 2.0, FS > 60.0"
  } else {
    // User option to choose which SNP filters will be used
    filter_ch = Channel.value("${params.VF_Filts}")
    println "User defined SNP filters set: ${params.VF_Filts}"
  }

  
  process GATK_VF {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.VF_Forks

    tag { chrom }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.VF_threads}:mem=${params.VF_memory}gb -lwalltime=${params.VF_walltime}:00:00"

    publishDir(
      path: "${params.VFDir}",
      mode: 'move',
    )

    input:
    //each chrom from chromosomes_ch
    tuple chrom, path(vcf), path(idx) from vf_ch
    val Filters from filter_ch
    path ref_genome
    path ref_dict
    path ref_index

    output:
    path("${chrom}.${params.VC_mode}.SNPs.filtered.vcf")
    path("${chrom}.${params.VC_mode}.SNPs.filtered.vcf.idx")
    
    beforeScript 'module load anaconda3/personal; source activate NF_GATK'

    script:
    """
    if [ ! -d tmp ]; then mkdir tmp; fi
    n_slots=`expr ${params.VF_threads} / 2 - 3`
    if [ \$n_slots -le 0 ]; then n_slots=1; fi
    taskset -c 0-\${n_slots} gatk VariantFiltration \\
      --tmp-dir tmp/ \\
      -R ${ref_genome} \\
      -V ${vcf} \\
      -O temp_filtered.vcf \\
      ${Filters} ${params.VF_args}

    taskset -c 0-\${n_slots} gatk SelectVariants \\
      --tmp-dir tmp/ \\
      -R ${ref_genome} \\
      -V temp_filtered.vcf \\
      -O ${chrom}.${params.VC_mode}.SNPs.filtered.vcf \\
      -select 'vc.isNotFiltered()'
    """
  }
}
