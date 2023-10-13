#!/bin/bash

# Estimate fraction of genome consisting of either variants types

# srun -p small -c 1 -J LD_SVs_SNPs -o log/LD_SVs_SNPs_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs.sh &

# parallel -a 02_infos/chr_list.txt -j 10 srun -p small -c 1 -J LD_SVs_SNPs_{} -o log/LD_SVs_SNPs_{}_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs.sh {} &

# VARIABLES
#CHR_LIST='02_infos/chr_list.txt'
CHR=$1

VCF_DIR="04_vcf"
SV_VCF="$VCF_DIR/SVs/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"
SNP_VCF="$VCF_DIR/SNPs/SNPs_MAF0.05_FMISS0.5.vcf.gz"

LD_DIR="LD/SVs_vs_SNPs"

MAX_DIST=100000

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load python/3.7


if [[ ! -d $LD_DIR ]]
then
  mkdir $LD_DIR
fi


# 1. Extract sites on given chromosome in both VCF files and output CHROM, POS, END, ID and GT fields to table
bcftools view -r $CHR $SV_VCF | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t[%GT\t]\n' > $VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_GTs_"$CHR".txt
bcftools view -r $CHR $SNP_VCF | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t[%GT\t]\n' > $VCF_DIR/SNPs/"$(basename -s .vcf.gz $SNP_VCF)"_GTs_"$CHR".txt

# 2. Run LD computation script on these files
python 01_scripts/utils/LD_SVs_SNPs.py $VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_GTs_"$CHR".txt $VCF_DIR/SNPs/"$(basename -s .vcf.gz $SNP_VCF)"_GTs_"$CHR".txt $MAX_DIST $LD_DIR/LD_SVs_SNPs_"$CHR".txt

# Remove files for each chromosome
#rm $VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_GTs_"$CHR".txt
#rm $VCF_DIR/SNPs/"$(basename -s .vcf.gz $SNP_VCF)"_GTs_"$CHR".txt
