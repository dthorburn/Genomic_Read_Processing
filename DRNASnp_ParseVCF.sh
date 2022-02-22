#!/bin/sh
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=4:mem=12gb
#PBS -N DRNASnp_ParseVCF
#PBS -j oe

module load anaconda3/personal
source activate NF_GATK

Project_Dir=/rds/general/user/dthorbur/home/ephemeral/05_RNA_GATKNF
cd $Project_Dir

## Automated handling of the difference in the DNAseq and RNAseq pipelines and variable numbers of files emitted from GATK_Variant_Call.nf
if [ -d ${Project_Dir}/07_MergedVCFs ]
then
	files=`ls -1 ${Project_Dir}/07_MergedVCFs/*.vcf`
	mode="RNAseq"
else
	files=`ls -1 ${Project_Dir}/06_SelectVariants/*SNPs.raw.vcf`
	mode="DNAseq"
fi
array=(${files})

## NB. If the VCFs are really big, this needs to be converted to an array. 
for vcf in "${array[@]}"
do
	echo $vcf
	Rscript ${Project_Dir}/DRNASnp_ParseVCF.R $vcf $mode
	echo "Done with ${vcf}"
done
