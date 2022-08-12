#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=8:mem=12gb
#PBS -N NF_GATK_Caller_Coordinator
#PBS -j oe

Project_Dir=/path/to/project/dir
cd $Project_Dir

module load nextflow/20.10.0

## To see help message paste the following command (without ##) into the console and press enter:
## module load nextflow/20.10.0; nextflow run GATK_Variant_Call.nf -c GATK_Variant_Call.config --profile imperial --help

echo "Starting: `date`"
nextflow run GATK_Variant_Call.nf \
	-c GATK_Variant_Call.config \
	--profile imperial \
	--VC_mode "DNAseq" \
	--RefGen "/path/to/genome.fasta" \
	--BamDir "/path/to/input/bams" 
echo "Finished: `date`"

