#!/bin/sh
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=4:mem=12gb
#PBS -N DRNASnp_ParseVCF
#PBS -j oe

module load anaconda3/personal
source activate NF_GATK

Project_Dir=/path/to/project/dir
cd $Project_Dir

ls -1 ${Project_Dir}/06_SelectVariants/*SNPs.raw.vcf > ${Project_Dir}/VCF_List.txt

## 1-7 as there are typically 7 chromosomes we are interested in for A. gambiae
## NB. If the VCFs are really big, this needs to be converted to an array. 
for i in {1..7..1}
do
	vcf=`sed ${i}"q;d" ${Project_Dir}/VCF_List.txt`
	Rscript ${Project_Dir}/DRNASnp_ParseVCF.R $vcf
	echo "Done with ${vcf}"
done
