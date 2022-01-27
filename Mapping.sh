#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=32:mem=62gb
#PBS -N NF_Mapping_Coordinator
#PBS -j oe

Project_Dir=/rds/general/user/dthorbur/ephemeral/02_Data/01_Pure_HFA
cd $Project_Dir

module load nextflow/20.10.0

echo "Starting: `date`"
nextflow run Mapping.nf -c nextflow.config --profile imperial --fastqc
echo "Finished: `date`"
