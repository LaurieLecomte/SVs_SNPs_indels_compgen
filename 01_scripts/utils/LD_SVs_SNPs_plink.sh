#!/bin/bash

# Estimate fraction of genome consisting of either variants types

# srun -p small -c 3 -J LD_SVs_SNPs_plink -o log/LD_SVs_SNPs_plink_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs_plink.sh &

# srun -p small -c 5 --mem=50G -J LD_SVs_SNPs_plink -o log/LD_SVs_SNPs_plink_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs_plink.sh &

# parallel -a 02_infos/chr_list.txt -j 10 srun -p small -c 3 -J LD_SVs_SNPs_plink_{} -o log/LD_SVs_SNPs_plink_{}_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs_plink.sh {} &

# VARIABLES
#CHR_LIST='02_infos/chr_list.txt'
CHR=$1

VCF_DIR="04_vcf"
SV_VCF="$VCF_DIR/SVs/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"
SNP_VCF="$VCF_DIR/SNPs/SNPs_MAF0.05_FMISS0.5.vcf.gz"


LD_DIR="LD/SVs_vs_SNPs"

MAX_DIST=100000

CPU=5

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load python/3.7
module load plink

if [[ ! -d $LD_DIR ]]
then
  mkdir $LD_DIR
fi


# 1. Extract sites on given chromosome in both VCF files and output CHROM, POS, END, ID and GT fields to table
#bcftools view -r $CHR $SV_VCF | bcftools sort -Oz > $VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_"$CHR".vcf.gz
#bcftools view -r $CHR $SNP_VCF | bcftools sort -Oz > $VCF_DIR/SNPs/"$(basename -s .vcf.gz $SNP_VCF)"_"$CHR".vcf.gz

# 2. Concatenate these 2 VCFs

#echo -e "$VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_"$CHR".vcf.gz\n$VCF_DIR/SNPs/"$(basename -s .vcf.gz $SNP_VCF)"_"$CHR".vcf.gz" > $LD_DIR/LD_SVs_SNPs_"$CHR"_list.txt

#bcftools concat -f $LD_DIR/LD_SVs_SNPs_"$CHR"_list.txt | bcftools sort -Oz > $LD_DIR/LD_SVs_SNPs_"$CHR".vcf.gz


# Run plink
#plink -r2 --allow-extra-chr --ld-window-kb 1000 --vcf $LD_DIR/LD_SVs_SNPs_"$CHR".vcf.gz --threads $CPU --out $LD_DIR/LD_SVs_SNPs_"$CHR"




# Replace POS by middle position between POS and END for SV larger than MAX_WIDTH
MAX_WIDTH=200
#Rscript 01_scripts/utils/recode_SVs_fromVCF.R $SV_VCF $VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_recodedPOS_"$MAX_WIDTH"bp.vcf $MAX_WIDTH

#bgzip $VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_recodedPOS_"$MAX_WIDTH"bp.vcf 
tabix -p vcf $VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_recodedPOS_"$MAX_WIDTH"bp.vcf.gz 

# Concatenate SVs and SNPs
echo -e "$VCF_DIR/SVs/"$(basename -s .vcf.gz $SV_VCF)"_recodedPOS_"$MAX_WIDTH"bp.vcf.gz\n$SNP_VCF" > $LD_DIR/LD_SVs_SNPs_list.txt

bcftools concat -f $LD_DIR/LD_SVs_SNPs_list.txt --threads $CPU -a | bcftools sort -Oz > $LD_DIR/SVs_SNPs_concat.vcf.gz


WIN_KB=10000

plink --r2 --allow-extra-chr --ld-window-kb $WIN_KB --ld-window-r2 0 --vcf $LD_DIR/SVs_SNPs_concat.vcf.gz --threads $CPU --out $LD_DIR/SVs_SNPs_concat_"$WIN_KB"









