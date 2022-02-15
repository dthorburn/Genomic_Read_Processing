#!/bin/sh
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=4:mem=12gb
#PBS -N DRNASnp_ParseVCF
#PBS -J 1-7
#PBS -j oe

module load anaconda3/personal
source activate NF_GATK

Project_Dir=/path/to/project/dir
cd $Project_Dir

ls -1 ${Project_Dir}/06_SelectVariants/*SNPs.raw.vcf > ${Project_Dir}/VCF_List.txt
vcf=`sed ${PBS_ARRAY_INDEX}"q;d" ${Project_Dir}/VCF_List.txt`

Rscript ${Project_Dir}/DRNASnp_ParseVCF.R $vcf
echo "Done with ${vcf}"
