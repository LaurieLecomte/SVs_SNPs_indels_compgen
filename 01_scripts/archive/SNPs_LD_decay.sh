#!/bin/bash 

# Get LD decay for SNPs

# Manitou
# srun -p small -c 1 --time=1:30:00 -J SNPs_LD_decay -o log/SNPs_LD_decay_%j.log /bin/sh 01_scripts/SNPs_LD_decay.sh  &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/SNPs"
RAW_SNPS_VCF="$RAW_VCF_DIR/SNPs_MAF0.05_FMISS0.5.vcf.gz"

ANGSD_INPUT_DIR="05_angsd_inputs/SNPs"
ANGSD_STATS_DIR="06_angsd_stats/SNPs"
ANGSD_BYPOP_DIR="07_angsd_bypop/SNPs"
ANGSD_FST_DIR="08_angsd_fst/SNPs"

PCA_DIR="09_pca/SNPs"
RDA_DIR="10_rda/SNPs"
#ANNOT_DIR="11_annotation/SNPs"
GO_DIR="12_go/SNPs"
FISHER_DIR="11_fisher_tests/SNPs"

ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SNPS_VCF_ANGSD=$RAW_SNPS_VCF

N_IND="$(less $ID_SEX_POP | wc -l)"

#MIN_MAF=0.05
#MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

LD_DIR="LD/SNPs"
CHR_LIST="02_infos/chr_list.txt"

# LOAD REQUIRED MODULES
module load PopLDecay/3.41
module load vcftools/0.1.16
module load bcftools/1.13

# 0. Create output dir
if [[ ! -d $LD_DIR ]]
then
  mkdir $LD_DIR
fi


# First split by chromosomes
bcftools index $SNPS_VCF_ANGSD

less $CHR_LIST | head -n2 | while read chr; 
do
  # Subset VCF by chromosomes
  bcftools filter -r $chr $SNPS_VCF_ANGSD -Oz > $LD_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)"."$chr".vcf.gz
  # Run PopLDdecay
  PopLDdecay -OutFilterSNP -InVCF $LD_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)"."$chr".vcf.gz -OutStat $LD_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)"."$chr".RO.stat -SubPop 02_infos/RO_IDs.txt -MAF $MIN_MAF -Miss $MAX_MISS -Het 1 -OutType 2
done
