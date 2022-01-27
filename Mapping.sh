#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=32:mem=62gb
#PBS -N NF_Mapping_Coordinator
#PBS -j oe

Project_Dir=/path/to/project/dir
cd $Project_Dir

module load nextflow/20.10.0

echo "Starting: `date`"
nextflow run Mapping.nf -c Mapping.config --profile imperial
echo "Finished: `date`"
