#!/bin/bash

# 

# srun -p small -c 3 -J LD_SVs_SNPs_plink -o log/LD_SVs_SNPs_plink_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs_plink.sh &

# srun -p small -c 3 --mem=200G -J LD_SVs_SNPs_plink -o log/LD_SVs_SNPs_plink_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs_plink.sh &

# parallel -a 02_infos/chr_list.txt -j 10 srun -p small -c 3 -J LD_SVs_SNPs_plink_{} -o log/LD_SVs_SNPs_plink_{}_%j.log /bin/sh 01_scripts/utils/LD_SVs_SNPs_plink.sh {} &

# VARIABLES
#CHR_LIST='02_infos/chr_list.txt'
CHR=$1

VCF_DIR="04_vcf"
SV_VCF="$VCF_DIR/SVs/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"
SNP_VCF="$VCF_DIR/SNPs/SNPs_MAF0.05_FMISS0.5.vcf.gz"


LD_DIR="LD"

MAX_DIST=100000

CPU=3

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load python/3.7
module load plink
module load bedtools

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


MAX_WIDTH=200


WIN_KB=100






# Try with no recoding of SV POS, and pre-remove SNPs overlapping with SVs

## Convert SNPs and SVs to bed
#

bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $SV_VCF > "${SV_VCF%.vcf.gz}".bed
#bcftools query -f '%CHROM\t%POS0\t%END\n' $SNP_VCF > "${SNP_VCF%.vcf.gz}".bed


## Remove SNPs overlapping with a SV
#bedtools window -a $SNP_VCF -b $SV_VCF -w 10 -v -header | bcftools sort -Oz > $LD_DIR/non_overlapping_SNPs_10bp.vcf.gz
#tabix -p vcf $LD_DIR/non_overlapping_SNPs_10bp.vcf.gz

## Concatenate non overlapping SNPs and SVs together
#echo -e "$SV_VCF\n$LD_DIR/non_overlapping_SNPs_10bp.vcf.gz" > $LD_DIR/LD_SVs_SNPs_list.txt
#bcftools concat -f $LD_DIR/LD_SVs_SNPs_list.txt --threads $CPU -a | bcftools sort -Oz > $LD_DIR/SVs_SNPs_concat_no_overlap.vcf.gz


#plink --r2 --allow-extra-chr --chr $CHR --ld-window 100 --ld-window-kb $WIN_KB --ld-window-r2 0 --vcf $LD_DIR/SVs_vs_SNPs/SVs_SNPs_concat_no_overlap.vcf.gz --threads $CPU --out $LD_DIR/SVs_SNPs_concat_"$WIN_KB"_$CHR


# For SNPs only 
plink --r2 --allow-extra-chr --chr $CHR --ld-window 100 --ld-window-kb $WIN_KB --ld-window-r2 0 --vcf $SNP_VCF --threads $CPU --out $LD_DIR/SNPs_vs_SNPs/SNPs_"$WIN_KB"_$CHR



