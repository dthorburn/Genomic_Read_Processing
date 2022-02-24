#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=8:mem=12gb
#PBS -N NF_Mapping_Coordinator
#PBS -j oe

Project_Dir=/path/to/project/dir
cd $Project_Dir

module load nextflow/20.10.0

echo "Starting: `date`"
nextflow run Mapping.nf -c Mapping.config \
	--profile imperial \
	--RefGen "/path/to/genome.{fasta,fa,fna}" \
	--InDir "/path/to/raw_fastq.gzs/" 
echo "Finished: `date`"
