#!/bin/bash

# Convert VCF to beagle and normalize genotype likelihoods

# manitou
# srun -p small -c 1 -J 01.3_SVs_vcf_to_beagle -o log/01.3_SVs_vcf_to_beagle_%j.log /bin/sh 01_scripts/01.3_SVs_vcf_to_beagle.sh &

# valeria
# srun -p ibis_small -c 1 -J 01.3_SVs_vcf_to_beagle -o log/01.3_SVs_vcf_to_beagle_%j.log /bin/sh 01_scripts/01.3_SVs_vcf_to_beagle.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/SVs"
RAW_SV_VCF="$RAW_VCF_DIR/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"

ANGSD_INPUT_DIR="05_angsd_inputs/SVs"
ANGSD_STATS_DIR="06_angsd_stats/SVs"
ANGSD_BYPOP_DIR="07_angsd_bypop/SVs"
ANGSD_FST_DIR="08_angsd_fst/SVs"

PCA_DIR="09_pca/SVs"
RDA_DIR="10_rda/SVs"
ANNOT_DIR="11_annotation/SVs"
GO_DIR="12_go/SVs"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SV_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=1

MIN_MAF=0.05
MAX_MAF=0.95

CHR_LIST="02_infos/chr_list.txt"

#$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz


# LOAD REQUIRED MODULES
module load bcftools/1.13
module load vcftools/0.1.16
module load R/4.1

 
# 1. Convert VCF to beagle, by looping on each chromosome
less $CHR_LIST | while read CHR
do 
  echo "convert vcf to beagle for $CHR"
  #vcftools --gzvcf $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz --BEAGLE-PL --chr $CHR --out $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF"_$CHR
  vcftools --gzvcf $SV_VCF_ANGSD --BEAGLE-PL --chr $CHR --out $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_"$CHR"
done

## Extract header from 1st chromosome beagle
FIRST_CHR="$(head -n1 $CHR_LIST)"
#head -n 1 $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF"_"$FIRST_CHR".BEAGLE.PL > $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".beagle
head -n 1 $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_"$FIRST_CHR".BEAGLE.PL > $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".beagle

## make a new file with the header
#cp "$VCF_ANGSD".beagleheader "$VCF_ANGSD".beagle

## remove header of individual beagles and append to global beagle
less $CHR_LIST | while read CHR
do 
  echo "append beagle from $CHR"
  #tail -n +2 $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF"_"$CHR".BEAGLE.PL >> $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".beagle
  tail -n +2 $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_"$CHR".BEAGLE.PL >> $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".beagle
done

## confirm SV count in beagle is = to number of SV in input VCF
#wc -l $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".beagle
wc -l $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".beagle


# 2. Normalize likelihoods in beagle file
#Rscript 01_scripts/utils/normalize_beagle.R $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".beagle $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle.contents
Rscript 01_scripts/utils/normalize_beagle.R $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".beagle $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".norm.beagle.contents


## add header
#head -n 1 $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF"_"$FIRST_CHR".BEAGLE.PL > $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle
head -n 1 $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_"$FIRST_CHR".BEAGLE.PL > $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".norm.beagle

#less $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle.contents >> $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle
less $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".norm.beagle.contents >> $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".norm.beagle

# 3. Compress 
#gzip $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle
gzip $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".norm.beagle

# Clean up
#for file in $(ls -1 $ANGSD_STATS_DIR/*.BEAGLE.PL); do rm $file; done
for file in $(ls -1 $ANGSD_INPUT_DIR/*.BEAGLE.PL); do rm $file; done
#for file in $(ls -1 $ANGSD_STATS_DIR/*.log); do rm $file; done
for file in $(ls -1 $ANGSD_INPUT_DIR/*.log); do rm $file; done
#rm $ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle.contents
rm $ANGSD_INPUT_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".norm.beagle.contents